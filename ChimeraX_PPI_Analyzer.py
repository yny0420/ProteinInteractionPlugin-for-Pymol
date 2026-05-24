"""
ChimeraX Protein Interaction Analyzer.

Load in ChimeraX:
    open /Users/yangyu/Desktop/interfacefinder/ChimeraX_PPI_Analyzer.py

Commands after loading:
    cxppi "model=#1 chains=A+B"
    cxppi_batch "#1:A+B;#2:A+B"
    cxppi_help

This is a ChimeraX command-script port of the PyMOL analyzer.  It uses a
Shrake-Rupley style SASA estimate to derive dASA/interface area and simple
distance-based rules for H-bonds, salt bridges, hydrophobic contacts, and
close contacts.
"""

from __future__ import print_function

import csv
import math
import os
import re


PLUGIN = "ChimeraX PPI Analyzer"

VDW_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80,
    "P": 1.80, "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98,
}

ACIDIC_RES = set(["ASP", "GLU"])
BASIC_RES = set(["ARG", "LYS", "HIS", "HID", "HIE", "HIP"])
HYDROPHOBIC_RES = set(["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TYR", "TRP", "PRO"])
POLAR_ELEMENTS = set(["N", "O", "S"])

ACIDIC_ATOMS = set(["OD1", "OD2", "OE1", "OE2"])
BASIC_ATOMS = set(["NZ", "NH1", "NH2", "NE", "ND1", "NE2"])

COLORS = {
    "chain_a": (80, 170, 255, 255),
    "chain_b": (240, 190, 110, 255),
    "interface_a": (0, 55, 180, 255),
    "interface_b": (255, 135, 0, 255),
    "hbond": (0, 255, 255, 255),
    "salt": (255, 0, 255, 255),
    "hydrophobic": (255, 230, 0, 255),
    "clash": (255, 0, 0, 255),
}


def _log(session, text):
    session.logger.info(text)


def _strip_quotes(text):
    text = str(text).strip()
    while len(text) >= 2 and ((text[0] == text[-1] == '"') or (text[0] == text[-1] == "'")):
        text = text[1:-1].strip()
    return text.strip('"').strip("'").strip()


def _safe_name(text, fallback="CXppi"):
    name = re.sub(r"[^A-Za-z0-9_]+", "_", str(text)).strip("_")
    if not name:
        name = fallback
    if not name[0].isalpha():
        name = "%s_%s" % (fallback, name)
    return name


def _element(atom):
    try:
        return atom.element.name.upper()
    except Exception:
        return str(getattr(atom, "element", "")).upper()


def _coord(atom):
    c = atom.coord
    return (float(c[0]), float(c[1]), float(c[2]))


def _dist(a, b):
    ca, cb = _coord(a), _coord(b)
    return math.sqrt((ca[0] - cb[0]) ** 2 + (ca[1] - cb[1]) ** 2 + (ca[2] - cb[2]) ** 2)


def _residue_name(atom):
    return atom.residue.name.upper()


def _chain_id(atom):
    return atom.residue.chain_id


def _resi(atom):
    return str(atom.residue.number)


def _residue_key(atom):
    return (_chain_id(atom), _resi(atom), _residue_name(atom))


def _residue_label(atom):
    return "%s:%s%s" % (_chain_id(atom), _residue_name(atom), _resi(atom))


def _atom_label(atom):
    return "%s.%s" % (_residue_label(atom), atom.name)


def _atom_radius(atom, probe=1.4):
    return VDW_RADII.get(_element(atom), 1.70) + float(probe)


def _heavy(atom):
    return _element(atom) != "H"


def _parse_key_values(spec):
    text = _strip_quotes(spec)
    out = {}
    for part in text.split():
        if "=" in part:
            key, value = part.split("=", 1)
            out[key.strip().lower()] = _strip_quotes(value)
    return out


def _parse_chains(chains):
    values = [x.strip() for x in re.split(r"[,+;:/\s]+", str(chains)) if x.strip()]
    if len(values) < 2:
        raise ValueError("Need at least two chains, e.g. chains=A+B")
    return values


def _parse_batch(spec):
    text = _strip_quotes(spec)
    jobs = []
    for item in [x.strip() for x in text.split(";") if x.strip()]:
        if ":" not in item:
            raise ValueError("Batch item must look like #1:A+B")
        model_spec, chains = item.split(":", 1)
        jobs.append((_strip_quotes(model_spec), _strip_quotes(chains)))
    return jobs


def _atomic_structures(session):
    from chimerax.atomic import AtomicStructure
    return list(session.models.list(type=AtomicStructure))


def _find_model(session, model_spec):
    text = _strip_quotes(model_spec)
    if text.startswith("#"):
        text = text[1:]
    models = _atomic_structures(session)
    for model in models:
        if getattr(model, "id_string", "") == text:
            return model
    for model in models:
        if getattr(model, "name", "") == text:
            return model
    raise ValueError("Cannot find atomic model: %s" % model_spec)


def _model_atoms(model):
    return [a for a in model.atoms if _heavy(a)]


def _chain_atoms(model, chain):
    chain = str(chain).strip()
    return [a for a in _model_atoms(model) if _chain_id(a) == chain]


def _fibonacci_sphere(n=64):
    pts = []
    offset = 2.0 / n
    inc = math.pi * (3.0 - math.sqrt(5.0))
    for i in range(n):
        y = ((i * offset) - 1.0) + (offset / 2.0)
        r = math.sqrt(max(0.0, 1.0 - y * y))
        phi = i * inc
        pts.append((math.cos(phi) * r, y, math.sin(phi) * r))
    return pts


SPHERE_POINTS = _fibonacci_sphere(72)


def _grid_key(coord, cell):
    return (int(math.floor(coord[0] / cell)), int(math.floor(coord[1] / cell)), int(math.floor(coord[2] / cell)))


def _build_grid(atoms, cell=6.0):
    grid = {}
    for atom in atoms:
        key = _grid_key(_coord(atom), cell)
        grid.setdefault(key, []).append(atom)
    return grid


def _nearby_atoms(atom, grid, cell=6.0):
    key = _grid_key(_coord(atom), cell)
    nearby = []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                nearby.extend(grid.get((key[0] + dx, key[1] + dy, key[2] + dz), []))
    return nearby


def sasa_by_atom(target_atoms, blocker_atoms, probe=1.4, points=None):
    points = points or SPHERE_POINTS
    grid = _build_grid(blocker_atoms)
    areas = {}
    for atom in target_atoms:
        center = _coord(atom)
        radius = _atom_radius(atom, probe)
        accessible = 0
        blockers = [b for b in _nearby_atoms(atom, grid) if b is not atom]
        for p in points:
            test = (center[0] + radius * p[0], center[1] + radius * p[1], center[2] + radius * p[2])
            buried = False
            for blocker in blockers:
                bc = _coord(blocker)
                br = _atom_radius(blocker, probe)
                if (test[0] - bc[0]) ** 2 + (test[1] - bc[1]) ** 2 + (test[2] - bc[2]) ** 2 < br ** 2:
                    buried = True
                    break
            if not buried:
                accessible += 1
        areas[atom] = 4.0 * math.pi * radius * radius * (float(accessible) / len(points))
    return areas


def calculate_dasa(model, chain_a, chain_b, cutoff=1.0):
    atoms_a = _chain_atoms(model, chain_a)
    atoms_b = _chain_atoms(model, chain_b)
    complex_atoms = atoms_a + atoms_b
    alone_a = sasa_by_atom(atoms_a, atoms_a)
    alone_b = sasa_by_atom(atoms_b, atoms_b)
    complex_area = sasa_by_atom(complex_atoms, complex_atoms)

    residues = {}
    for atom in complex_atoms:
        alone = alone_a.get(atom, alone_b.get(atom, 0.0))
        diff = max(0.0, alone - complex_area.get(atom, 0.0))
        key = _residue_key(atom)
        row = residues.setdefault(key, {
            "chain": key[0], "resi": key[1], "resn": key[2],
            "dasa": 0.0, "atoms": [],
        })
        row["dasa"] += diff
        row["atoms"].append(atom)

    selected = [r for r in residues.values() if r["dasa"] >= float(cutoff)]
    selected.sort(key=lambda r: r["dasa"], reverse=True)
    total_buried = sum(r["dasa"] for r in selected)
    return {
        "residues": selected,
        "total_buried_area": total_buried,
        "interface_area": total_buried / 2.0,
    }


def _closest_pairs(atoms_a, atoms_b, cutoff, by_residue=False, predicate=None):
    hits = {}
    for a in atoms_a:
        for b in atoms_b:
            if predicate and not predicate(a, b):
                continue
            d = _dist(a, b)
            if d <= cutoff:
                key = (_residue_key(a), _residue_key(b)) if by_residue else (a, b)
                old = hits.get(key)
                if old is None or d < old["distance"]:
                    hits[key] = {"atom_a": a, "atom_b": b, "distance": d}
    return sorted(hits.values(), key=lambda x: x["distance"])


def _is_polar_pair(a, b):
    return _element(a) in POLAR_ELEMENTS and _element(b) in POLAR_ELEMENTS


def _is_salt_pair(a, b):
    acidic_a = _residue_name(a) in ACIDIC_RES and a.name in ACIDIC_ATOMS
    basic_a = _residue_name(a) in BASIC_RES and a.name in BASIC_ATOMS
    acidic_b = _residue_name(b) in ACIDIC_RES and b.name in ACIDIC_ATOMS
    basic_b = _residue_name(b) in BASIC_RES and b.name in BASIC_ATOMS
    return (acidic_a and basic_b) or (basic_a and acidic_b)


def _is_hydrophobic_pair(a, b):
    return (
        _residue_name(a) in HYDROPHOBIC_RES and _residue_name(b) in HYDROPHOBIC_RES
        and _element(a) in ("C", "S") and _element(b) in ("C", "S")
        and a.name not in ("C", "CA") and b.name not in ("C", "CA")
    )


def find_contacts(model, chain_a, chain_b, hbond_cutoff=3.6, salt_cutoff=4.0,
                  hydrophobic_cutoff=4.2, clash_cutoff=2.2):
    atoms_a = _chain_atoms(model, chain_a)
    atoms_b = _chain_atoms(model, chain_b)
    return {
        "hbonds": _closest_pairs(atoms_a, atoms_b, hbond_cutoff, False, _is_polar_pair),
        "salt_bridges": _closest_pairs(atoms_a, atoms_b, salt_cutoff, True, _is_salt_pair),
        "hydrophobic": _closest_pairs(atoms_a, atoms_b, hydrophobic_cutoff, True, _is_hydrophobic_pair),
        "clashes": _closest_pairs(atoms_a, atoms_b, clash_cutoff, False, None),
    }


def score_breakdown(interface_area, contacts):
    area = min(float(interface_area) / 1200.0 * 40.0, 40.0)
    hbonds = min(len(contacts["hbonds"]) * 2.0, 20.0)
    salt = min(len(contacts["salt_bridges"]) * 4.0, 16.0)
    hydrophobic = min(len(contacts["hydrophobic"]) * 0.35, 20.0)
    clash_penalty = min(len(contacts["clashes"]) * 3.0, 30.0)
    total = max(0.0, min(100.0, area + hbonds + salt + hydrophobic - clash_penalty))
    return {
        "area": area, "hbonds": hbonds, "salt": salt,
        "hydrophobic": hydrophobic, "clash_penalty": clash_penalty,
        "total": total,
    }


def local_area_rows(interface):
    rows = []
    total = interface["total_buried_area"]
    for r in interface["residues"]:
        rows.append({
            "chain": r["chain"], "residue": "%s%s" % (r["resn"], r["resi"]),
            "dasa": r["dasa"], "local_area": r["dasa"] / 2.0,
            "share": r["dasa"] / total if total else 0.0,
        })
    return rows


def color_atoms(atoms, rgba):
    for a in atoms:
        try:
            a.color = rgba
        except Exception:
            pass


def show_atoms(atoms):
    for a in atoms:
        try:
            a.display = True
            a.residue.ribbon_display = True
        except Exception:
            pass


def _pseudobond_group(session, name, color):
    group = session.pb_manager.get_group(name)
    group.clear()
    group.color = color
    group.radius = 0.06
    return group


def add_pseudobonds(session, name, contacts, color, max_items=120):
    group = _pseudobond_group(session, name, color)
    for c in contacts[:max_items]:
        pb = group.new_pseudobond(c["atom_a"], c["atom_b"])
        pb.length = c["distance"]
        pb.color = color
    return group


def visualize(session, model, chain_a, chain_b, interface, contacts, prefix):
    atoms_a = _chain_atoms(model, chain_a)
    atoms_b = _chain_atoms(model, chain_b)
    show_atoms(atoms_a + atoms_b)
    color_atoms(atoms_a, COLORS["chain_a"])
    color_atoms(atoms_b, COLORS["chain_b"])
    for r in interface["residues"]:
        color_atoms(r["atoms"], COLORS["interface_a"] if r["chain"] == chain_a else COLORS["interface_b"])
        show_atoms(r["atoms"])
    add_pseudobonds(session, prefix + "_hbonds", contacts["hbonds"], COLORS["hbond"])
    add_pseudobonds(session, prefix + "_salt", contacts["salt_bridges"], COLORS["salt"])
    add_pseudobonds(session, prefix + "_hydrophobic", contacts["hydrophobic"], COLORS["hydrophobic"], 80)
    add_pseudobonds(session, prefix + "_clashes", contacts["clashes"], COLORS["clash"], 80)


def contact_rows(kind, contacts, interface):
    lookup = {(r["chain"], "%s%s" % (r["resn"], r["resi"])): r["dasa"] for r in interface["residues"]}
    rows = []
    for c in contacts:
        a, b = c["atom_a"], c["atom_b"]
        ra, rb = _residue_label(a).split(":", 1), _residue_label(b).split(":", 1)
        da = lookup.get((ra[0], ra[1]), 0.0)
        db = lookup.get((rb[0], rb[1]), 0.0)
        rows.append({
            "type": kind,
            "residue_a": _residue_label(a),
            "residue_b": _residue_label(b),
            "atom_a": _atom_label(a),
            "atom_b": _atom_label(b),
            "distance": c["distance"],
            "local_pair_area": (da + db) / 2.0,
        })
    return rows


def build_table(model, chain_a, chain_b, interface, contacts, score):
    rows = []
    rows += contact_rows("H-bond", contacts["hbonds"], interface)
    rows += contact_rows("Salt bridge", contacts["salt_bridges"], interface)
    rows += contact_rows("Hydrophobic", contacts["hydrophobic"], interface)
    rows += contact_rows("Close contact", contacts["clashes"], interface)
    lines = []
    lines.append("[%s] %s chain %s vs chain %s" % (PLUGIN, model.name, chain_a, chain_b))
    lines.append("Whole interaction area: %.1f A^2 | Total buried area: %.1f A^2 | Score: %.1f/100" %
                 (interface["interface_area"], interface["total_buried_area"], score["total"]))
    lines.append("Score components: area %.1f + H-bond %.1f + salt %.1f + hydrophobic %.1f - clash %.1f" %
                 (score["area"], score["hbonds"], score["salt"], score["hydrophobic"], score["clash_penalty"]))
    lines.append("H-bonds: %d | Salt bridges: %d | Hydrophobic: %d | Close contacts: %d" %
                 (len(contacts["hbonds"]), len(contacts["salt_bridges"]), len(contacts["hydrophobic"]), len(contacts["clashes"])))
    lines.append("")
    lines.append("%-14s %-12s %-12s %8s %12s  %s" %
                 ("Type", "Residue A", "Residue B", "Dist(A)", "LocalArea", "Atoms"))
    lines.append("-" * 96)
    for r in rows[:250]:
        lines.append("%-14s %-12s %-12s %8.2f %12.1f  %s -- %s" %
                     (r["type"], r["residue_a"], r["residue_b"], r["distance"],
                      r["local_pair_area"], r["atom_a"], r["atom_b"]))
    lines.append("")
    lines.append("Local interface area hot spots")
    lines.append("%-8s %-10s %12s %14s %8s" % ("Chain", "Residue", "dASA(A2)", "LocalArea", "Share"))
    for r in local_area_rows(interface)[:30]:
        lines.append("%-8s %-10s %12.1f %14.1f %7.1f%%" %
                     (r["chain"], r["residue"], r["dasa"], r["local_area"], r["share"] * 100.0))
    return "\n".join(lines), rows


def analyze_one(session, model_spec, chains, cutoff=1.0, visualize_output=True,
                hbond_cutoff=3.6, salt_cutoff=4.0, hydrophobic_cutoff=4.2,
                clash_cutoff=2.2, prefix="CXppi"):
    model = _find_model(session, model_spec)
    chain_a, chain_b = _parse_chains(chains)[:2]
    interface = calculate_dasa(model, chain_a, chain_b, cutoff)
    contacts = find_contacts(model, chain_a, chain_b, hbond_cutoff, salt_cutoff, hydrophobic_cutoff, clash_cutoff)
    score = score_breakdown(interface["interface_area"], contacts)
    pfx = _safe_name("%s%s%s%s" % (prefix, getattr(model, "id_string", model.name), chain_a, chain_b))
    if visualize_output:
        visualize(session, model, chain_a, chain_b, interface, contacts, pfx)
    table, rows = build_table(model, chain_a, chain_b, interface, contacts, score)
    _log(session, table)
    return {
        "model": model, "chain_a": chain_a, "chain_b": chain_b,
        "interface": interface, "contacts": contacts, "score": score,
        "table": table, "rows": rows, "prefix": pfx,
    }


def cxppi(session, spec):
    opts = _parse_key_values(spec)
    model = opts.get("model") or opts.get("object") or opts.get("m")
    chains = opts.get("chains") or opts.get("chain") or opts.get("c")
    if not model or not chains:
        raise ValueError('Use: cxppi "model=#1 chains=A+B"')
    return analyze_one(
        session, model, chains,
        cutoff=float(opts.get("cutoff", 1.0)),
        visualize_output=opts.get("visualize", "true").lower() not in ("0", "false", "no"),
        hbond_cutoff=float(opts.get("hbond_cutoff", 3.6)),
        salt_cutoff=float(opts.get("salt_cutoff", 4.0)),
        hydrophobic_cutoff=float(opts.get("hydrophobic_cutoff", 4.2)),
        clash_cutoff=float(opts.get("clash_cutoff", 2.2)),
        prefix=opts.get("prefix", "CXppi"),
    )


def cxppi_batch(session, spec):
    jobs = _parse_batch(spec)
    lines = []
    rows = []
    lines.append("[%s] Batch summary" % PLUGIN)
    lines.append("%-10s %-8s %10s %8s %6s %6s %12s %8s %7s  %s" %
                 ("Model", "Pair", "Area(A2)", "Resid", "HB", "Salt", "Hydrophobic", "Clash", "Score", "Top local residues"))
    lines.append("-" * 120)
    for model_spec, chains in jobs:
        result = analyze_one(session, model_spec, chains, visualize_output=True)
        hotspots = "; ".join("%s:%s %.1fA2" % (r["chain"], r["residue"], r["local_area"])
                             for r in local_area_rows(result["interface"])[:5])
        row = {
            "model": result["model"].name,
            "pair": "%s-%s" % (result["chain_a"], result["chain_b"]),
            "area": result["interface"]["interface_area"],
            "residues": len(result["interface"]["residues"]),
            "hbonds": len(result["contacts"]["hbonds"]),
            "salt": len(result["contacts"]["salt_bridges"]),
            "hydrophobic": len(result["contacts"]["hydrophobic"]),
            "clashes": len(result["contacts"]["clashes"]),
            "score": result["score"]["total"],
            "hotspots": hotspots,
        }
        rows.append(row)
        lines.append("%-10s %-8s %10.1f %8d %6d %6d %12d %8d %7.1f  %s" %
                     (row["model"], row["pair"], row["area"], row["residues"],
                      row["hbonds"], row["salt"], row["hydrophobic"], row["clashes"],
                      row["score"], row["hotspots"]))
    table = "\n".join(lines)
    _log(session, table)
    return rows


def cxppi_help(session):
    text = """
ChimeraX PPI Analyzer commands:

  cxppi "model=#1 chains=A+B"
  cxppi "model=#1 chains=A+B cutoff=1.0 prefix=test"
  cxppi_batch "#1:A+B;#2:A+B;#3:C+D"

Output:
  Whole interaction area = total buried area / 2
  LocalArea per residue = residue dASA / 2
  Score = area + H-bond + salt + hydrophobic - clash penalty

This is a geometry-based heuristic, not a binding free-energy calculation.
"""
    _log(session, text)
    return text


def register_commands(session):
    from chimerax.core.commands import CmdDesc, RestOfLine, register
    register("cxppi", CmdDesc(required=[("spec", RestOfLine)]), cxppi, logger=session.logger)
    register("cxppi_batch", CmdDesc(required=[("spec", RestOfLine)]), cxppi_batch, logger=session.logger)
    register("cxppi_help", CmdDesc(), cxppi_help, logger=session.logger)
    _log(session, "%s loaded. Run: cxppi_help" % PLUGIN)


if "session" in globals():
    register_commands(session)

