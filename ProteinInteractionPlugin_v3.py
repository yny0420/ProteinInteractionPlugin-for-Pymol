"""
Protein Interaction Analyzer for PyMOL.

Copyright (c) yangyu. All rights reserved.

Load in PyMOL:
    run /Users/yangyu/Desktop/interfacefinder/ProteinInteractionPlugin.py

Commands:
    ppi_panel
    ppi_analyze 1brs, A+B
    ppi_analyze 1brs, A+B+C+D
    ppi_batch "obj1:A+B;obj2:C+D"
    ppi_fetch 1brs

The dASA interface calculation is adapted from the user's
InterfaceResidues_new.py workflow.
"""

from __future__ import print_function

import csv
import itertools
import math
import os
import re

from pymol import cmd, stored


PLUGIN = "Protein Interaction Analyzer"

ACIDIC_ATOMS = "(resn ASP+GLU and name OD1+OD2+OE1+OE2)"
BASIC_ATOMS = "(resn ARG+LYS+HIS+HID+HIE+HIP and name NZ+NH1+NH2+NE+ND1+NE2)"
POLAR_ATOMS = "(donors or acceptors)"
HYDROPHOBIC_ATOMS = "(resn ALA+VAL+LEU+ILE+MET+PHE+TYR+TRP+PRO and elem C+S and not name C+CA)"
HEAVY_ATOMS = "(not hydro)"

COLORS = {
    "chain_a": "lightblue",
    "chain_b": "wheat",
    "interface_a": "marine",
    "interface_b": "orange",
    "hbond": "cyan",
    "salt": "magenta",
    "hydrophobic": "yellow",
    "clash": "red",
}


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "HID": "H",
    "HIE": "H", "HIP": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "MSE": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def _aa1(resn):
    return AA3_TO_1.get(str(resn).upper(), str(resn)[:1])


def _safe_name(text, fallback="PPI"):
    name = re.sub(r"[^A-Za-z0-9]+", "", str(text))
    if not name:
        name = fallback
    if not name[0].isalpha():
        name = fallback + name
    return name


def _strip_quotes(text):
    text = str(text).strip()
    while len(text) >= 2 and (
        (text[0] == text[-1] == '"') or
        (text[0] == text[-1] == "'")
    ):
        text = text[1:-1].strip()
    return text.strip('"').strip("'").strip()


def _as_bool(value):
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() not in ("0", "false", "no", "off", "")


def _model_selection(obj_or_sel):
    text = _strip_quotes(obj_or_sel)
    if not text:
        raise ValueError("Object/selection cannot be empty")
    if text.lower() == "all":
        return "all"
    try:
        if text in cmd.get_names("objects"):
            return "(model %s)" % text
    except Exception:
        pass
    if text.startswith("(") and text.endswith(")"):
        return text
    return "(%s)" % text


def _chain_selection(chain):
    text = str(chain).strip()
    if not text:
        raise ValueError("Chain cannot be empty")
    low = text.lower()
    if low.startswith(("chain ", "c.", "chain.")):
        return text
    if text in ("*", "all"):
        return "all"
    return "chain %s" % text


def _parse_chains(chains):
    if isinstance(chains, (list, tuple)):
        raw = [str(c).strip() for c in chains]
    else:
        raw = [c.strip() for c in re.split(r"[,+;:\s]+", str(chains)) if c.strip()]
    values, seen = [], set()
    for chain in raw:
        if chain not in seen:
            values.append(chain)
            seen.add(chain)
    if len(values) < 2:
        raise ValueError("Please provide at least two chains, for example A+B")
    return values


def _parse_inline_command(first_arg, default_chains):
    """Handle PyMOL builds that pass 'object, chains, key=value' as one string."""
    text = _strip_quotes(first_arg)
    if "," not in text:
        return first_arg, default_chains, {}
    parts = [part.strip() for part in text.split(",") if part.strip()]
    if not parts:
        return first_arg, default_chains, {}
    complex_obj = _strip_quotes(parts[0])
    chains = _strip_quotes(parts[1]) if len(parts) > 1 and "=" not in parts[1] else default_chains
    options = {}
    for part in parts[1:]:
        if "=" not in part:
            continue
        key, value = [_strip_quotes(x) for x in part.split("=", 1)]
        options[key] = value
    return complex_obj, chains, options


def _parse_batch_specs(specs):
    """Parse 'obj1:A+B;obj2:C+D' or 'obj1:A+B|obj2:C+D'."""
    text = _strip_quotes(specs)
    if not text:
        raise ValueError("Please provide batch specs, e.g. 001:A+B;002:A+B")
    items = [x.strip() for x in re.split(r"[;|]+", text) if x.strip()]
    parsed = []
    for item in items:
        if ":" in item:
            obj, chains = [_strip_quotes(x) for x in item.split(":", 1)]
        elif "," in item:
            obj, chains = [_strip_quotes(x) for x in item.split(",", 1)]
        else:
            raise ValueError("Batch item must look like object:A+B, got: %s" % item)
        parsed.append((obj, chains))
    return parsed


def _float_option(options, key, current):
    if key not in options:
        return current
    return float(options[key])


def _string_option(options, key, current):
    return options.get(key, current)


def _pair_prefix(prefix, chain_a, chain_b):
    return _safe_name("%s%s%s" % (prefix, chain_a, chain_b), "PPI")


def _selref(name):
    return "(%%%s)" % name


def _dist(coord_a, coord_b):
    return math.sqrt(
        (coord_a[0] - coord_b[0]) ** 2 +
        (coord_a[1] - coord_b[1]) ** 2 +
        (coord_a[2] - coord_b[2]) ** 2
    )


def _atom_sel(atom):
    model = getattr(atom, "model", None) or getattr(atom, "_ppi_model", None)
    if model:
        return "(model %s and index %d)" % (model, atom.index)
    return "(index %d)" % atom.index


def _atom_label(atom):
    return "%s:%s%s.%s" % (atom.chain or "-", atom.resn, atom.resi, atom.name)


def _residue_label(atom):
    return "%s:%s%s" % (atom.chain or "-", atom.resn, atom.resi)


def _selection_model_name(selection):
    text = str(selection).strip()
    match = re.search(r"model\s+([^\s\)]+)", text)
    if match:
        return match.group(1)
    names = cmd.get_names("objects")
    for name in names:
        if re.search(r"(^|[^A-Za-z0-9_])%s([^A-Za-z0-9_]|$)" % re.escape(name), text):
            return name
    return None


def _atoms(selection):
    model_name = _selection_model_name(selection)
    atoms = list(cmd.get_model(selection).atom)
    for atom in atoms:
        if not hasattr(atom, "model") and model_name:
            atom._ppi_model = model_name
    return atoms


def _safe_delete(names):
    for name in names:
        try:
            cmd.delete(name)
        except Exception:
            pass


def _resi_sort_key(resi):
    digits = "".join(ch for ch in str(resi) if ch.isdigit() or ch == "-")
    suffix = "".join(ch for ch in str(resi) if not (ch.isdigit() or ch == "-"))
    try:
        return (int(digits), suffix)
    except Exception:
        return (999999, str(resi))


def ppi_fetch(pdb_id):
    """Fetch a PDB by ID, keeping object naming predictable."""
    pdb_id = str(pdb_id).strip()
    if not pdb_id:
        raise ValueError("Please provide a PDB ID")
    cmd.fetch(pdb_id.lower(), async_=0)
    return pdb_id.lower()


def ppi_load(path, object_name=""):
    """Load a local PDB/mmCIF file."""
    expanded = os.path.expanduser(str(path).strip())
    if not os.path.exists(expanded):
        raise IOError("File does not exist: %s" % expanded)
    name = _safe_name(object_name or os.path.splitext(os.path.basename(expanded))[0], "PPIObj")
    cmd.load(expanded, name)
    return name


def calculate_interface_dasa(complex_obj, chain_a, chain_b, cutoff=1.0, prefix="PPI"):
    """Create an interface selection and return residue dASA data."""
    complex_sel = _model_selection(complex_obj)
    sel_a = _chain_selection(chain_a)
    sel_b = _chain_selection(chain_b)
    cutoff = float(cutoff)

    old_dot = cmd.get("dot_solvent")
    old_solvent = cmd.get("solvent_radius")
    tmp_complex = cmd.get_unused_name("PPITmpComplex")
    tmp_a = cmd.get_unused_name("PPITmpA")
    tmp_b = cmd.get_unused_name("PPITmpB")
    tmp_sel = cmd.get_unused_name("PPITmpSel")
    interface_sel = _safe_name(prefix + "iface")
    residues = {}

    try:
        cmd.set("dot_solvent", 1)
        cmd.set("solvent_radius", 1.4)

        cmd.create(tmp_complex, complex_sel)
        cmd.remove("%s and not (polymer.protein and ((%s) or (%s)))" %
                   (tmp_complex, sel_a, sel_b))

        cmd.get_area(tmp_complex, load_b=1)
        cmd.alter(tmp_complex, "q=b")

        cmd.create(tmp_a, "%s and (%s)" % (tmp_complex, sel_a))
        cmd.create(tmp_b, "%s and (%s)" % (tmp_complex, sel_b))
        cmd.get_area(tmp_a, load_b=1)
        cmd.get_area(tmp_b, load_b=1)
        cmd.alter("(%s or %s)" % (tmp_a, tmp_b), "b=b-q")

        stored.ppi_dasa_atoms = []
        cmd.iterate(
            "(%s or %s)" % (tmp_a, tmp_b),
            "stored.ppi_dasa_atoms.append((chain, resi, resn, b))"
        )

        for chain, resi, resn, atom_dasa in stored.ppi_dasa_atoms:
            key = (chain, resi, resn)
            row = residues.setdefault(key, {
                "chain": chain,
                "resi": resi,
                "resn": resn,
                "dasa": 0.0,
                "positive_dasa": 0.0,
            })
            row["dasa"] += float(atom_dasa)
            if atom_dasa > 0:
                row["positive_dasa"] += float(atom_dasa)

        selected = [r for r in residues.values() if r["positive_dasa"] >= cutoff]
        selected.sort(key=lambda r: (r["chain"], _resi_sort_key(r["resi"])))

        cmd.select(tmp_sel, "none")
        for residue in selected:
            cmd.select(
                tmp_sel,
                "%s or (%s and chain %s and resi %s)" %
                (_selref(tmp_sel), complex_sel, residue["chain"], residue["resi"])
            )

        cmd.select(interface_sel, _selref(tmp_sel))
        cmd.enable(interface_sel)

        total_buried = sum(r["positive_dasa"] for r in selected)
        return {
            "selection": interface_sel,
            "residues": selected,
            "total_buried_area": total_buried,
            "interface_area": total_buried / 2.0,
        }
    finally:
        cmd.set("dot_solvent", old_dot)
        cmd.set("solvent_radius", old_solvent)
        stored.ppi_dasa_atoms = []
        _safe_delete([tmp_complex, tmp_a, tmp_b, tmp_sel])


def _closest_pairs(atoms_a, atoms_b, cutoff, by_residue=False):
    cutoff = float(cutoff)
    hits = {}
    for atom_a in atoms_a:
        for atom_b in atoms_b:
            distance = _dist(atom_a.coord, atom_b.coord)
            if distance <= cutoff:
                if by_residue:
                    key = (atom_a.chain, atom_a.resi, atom_a.resn,
                           atom_b.chain, atom_b.resi, atom_b.resn)
                else:
                    key = (atom_a.model, atom_a.index, atom_b.model, atom_b.index)
                old = hits.get(key)
                if old is None or distance < old["distance"]:
                    hits[key] = {
                        "atom_a": atom_a,
                        "atom_b": atom_b,
                        "distance": distance,
                    }
    return sorted(hits.values(), key=lambda x: x["distance"])


def find_hbonds(complex_obj, chain_a, chain_b, cutoff=3.6, angle=45):
    complex_sel = _model_selection(complex_obj)
    sel_a = "(%s and (%s) and %s)" % (complex_sel, _chain_selection(chain_a), POLAR_ATOMS)
    sel_b = "(%s and (%s) and %s)" % (complex_sel, _chain_selection(chain_b), POLAR_ATOMS)
    pairs = cmd.find_pairs(sel_a, sel_b, mode=1, cutoff=float(cutoff), angle=float(angle))
    contacts = []
    for left, right in pairs:
        model_a, index_a = left
        model_b, index_b = right
        atoms_a = _atoms("(model %s and index %d)" % (model_a, index_a))
        atoms_b = _atoms("(model %s and index %d)" % (model_b, index_b))
        if atoms_a and atoms_b:
            contacts.append({
                "atom_a": atoms_a[0],
                "atom_b": atoms_b[0],
                "distance": _dist(atoms_a[0].coord, atoms_b[0].coord),
            })
    return sorted(contacts, key=lambda x: x["distance"])


def find_salt_bridges(complex_obj, chain_a, chain_b, cutoff=4.0):
    complex_sel = _model_selection(complex_obj)
    acid_a = _atoms("(%s and (%s) and %s)" % (complex_sel, _chain_selection(chain_a), ACIDIC_ATOMS))
    base_a = _atoms("(%s and (%s) and %s)" % (complex_sel, _chain_selection(chain_a), BASIC_ATOMS))
    acid_b = _atoms("(%s and (%s) and %s)" % (complex_sel, _chain_selection(chain_b), ACIDIC_ATOMS))
    base_b = _atoms("(%s and (%s) and %s)" % (complex_sel, _chain_selection(chain_b), BASIC_ATOMS))
    contacts = _closest_pairs(acid_a, base_b, cutoff, True)
    contacts.extend(_closest_pairs(base_a, acid_b, cutoff, True))
    return sorted(contacts, key=lambda x: x["distance"])


def find_hydrophobic_contacts(complex_obj, chain_a, chain_b, cutoff=4.2):
    complex_sel = _model_selection(complex_obj)
    atoms_a = _atoms("(%s and (%s) and %s)" % (complex_sel, _chain_selection(chain_a), HYDROPHOBIC_ATOMS))
    atoms_b = _atoms("(%s and (%s) and %s)" % (complex_sel, _chain_selection(chain_b), HYDROPHOBIC_ATOMS))
    return _closest_pairs(atoms_a, atoms_b, cutoff, True)


def find_close_contacts(complex_obj, chain_a, chain_b, cutoff=2.2):
    complex_sel = _model_selection(complex_obj)
    atoms_a = _atoms("(%s and (%s) and %s)" % (complex_sel, _chain_selection(chain_a), HEAVY_ATOMS))
    atoms_b = _atoms("(%s and (%s) and %s)" % (complex_sel, _chain_selection(chain_b), HEAVY_ATOMS))
    return _closest_pairs(atoms_a, atoms_b, cutoff, False)


def score_interface(interface_area, contacts):
    """Heuristic score for comparing related structures, not binding energy."""
    score = min(float(interface_area) / 1200.0 * 40.0, 40.0)
    score += min(len(contacts["hbonds"]) * 2.0, 20.0)
    score += min(len(contacts["salt_bridges"]) * 4.0, 16.0)
    score += min(len(contacts["hydrophobic"]) * 0.35, 20.0)
    score -= min(len(contacts["clashes"]) * 3.0, 30.0)
    return max(0.0, min(100.0, score))


def score_breakdown(interface_area, contacts):
    """Return additive score components so interfaces can be compared transparently."""
    area = min(float(interface_area) / 1200.0 * 40.0, 40.0)
    hbonds = min(len(contacts["hbonds"]) * 2.0, 20.0)
    salt = min(len(contacts["salt_bridges"]) * 4.0, 16.0)
    hydrophobic = min(len(contacts["hydrophobic"]) * 0.35, 20.0)
    clash_penalty = min(len(contacts["clashes"]) * 3.0, 30.0)
    total = max(0.0, min(100.0, area + hbonds + salt + hydrophobic - clash_penalty))
    return {
        "area": area,
        "hbonds": hbonds,
        "salt_bridges": salt,
        "hydrophobic": hydrophobic,
        "clash_penalty": clash_penalty,
        "total": total,
    }


def local_area_rows(interface):
    """Residue-level area contributions derived from dASA."""
    rows = []
    for residue in interface["residues"]:
        dasa = residue["positive_dasa"]
        rows.append({
            "chain": residue["chain"],
            "residue": "%s%s" % (residue["resn"], residue["resi"]),
            "dasa": dasa,
            "local_interface_area": dasa / 2.0,
            "interface_fraction": (dasa / interface["total_buried_area"])
            if interface["total_buried_area"] else 0.0,
        })
    return sorted(rows, key=lambda r: r["dasa"], reverse=True)


def residue_area_lookup(interface):
    lookup = {}
    for residue in interface["residues"]:
        lookup[(residue["chain"], "%s%s" % (residue["resn"], residue["resi"]))] = residue["positive_dasa"]
    return lookup


def _distance_object(name, contacts, color, max_items=120):
    obj_name = _safe_name(name)
    cmd.delete(obj_name)
    if not contacts:
        return None
    for contact in contacts[:int(max_items)]:
        cmd.distance(obj_name, _atom_sel(contact["atom_a"]), _atom_sel(contact["atom_b"]))
    cmd.color(color, obj_name)
    cmd.set("dash_radius", 0.08, obj_name)
    cmd.set("dash_gap", 0.25, obj_name)
    cmd.set("label_size", 14, obj_name)
    cmd.set("label_color", color, obj_name)
    return obj_name


def make_visual_objects(complex_obj, chain_a, chain_b, interface, contacts, prefix):
    complex_sel = _model_selection(complex_obj)
    chain_a_sel = _safe_name(prefix + "chain" + str(chain_a))
    chain_b_sel = _safe_name(prefix + "chain" + str(chain_b))

    cmd.select(chain_a_sel, "%s and (%s)" % (complex_sel, _chain_selection(chain_a)))
    cmd.select(chain_b_sel, "%s and (%s)" % (complex_sel, _chain_selection(chain_b)))

    cmd.show("cartoon", _selref(chain_a_sel))
    cmd.show("cartoon", _selref(chain_b_sel))
    cmd.hide("lines", _selref(chain_a_sel))
    cmd.hide("lines", _selref(chain_b_sel))
    cmd.color(COLORS["chain_a"], _selref(chain_a_sel))
    cmd.color(COLORS["chain_b"], _selref(chain_b_sel))

    iface_sel = interface["selection"]
    iface_ref = _selref(iface_sel)
    cmd.show("sticks", iface_ref)
    cmd.color(COLORS["interface_a"], "%s and (%s)" % (iface_ref, _chain_selection(chain_a)))
    cmd.color(COLORS["interface_b"], "%s and (%s)" % (iface_ref, _chain_selection(chain_b)))
    cmd.util.cnc(iface_ref)
    residue_label_sel = "%s and name CA" % iface_ref
    cmd.label(residue_label_sel, '"%s:%s%s" % (chain, _aa1(resn), resi)')
    cmd.show("labels", residue_label_sel)
    cmd.set("label_color", "white", iface_sel)
    cmd.set("label_size", 18, iface_sel)
    cmd.set("label_outline_color", "black", iface_sel)
    cmd.set("label_position", [0.0, 0.0, 1.6], iface_sel)

    _distance_object(prefix + "hbonds", contacts["hbonds"], COLORS["hbond"], 160)
    _distance_object(prefix + "salt", contacts["salt_bridges"], COLORS["salt"], 160)
    _distance_object(prefix + "hydrophobic", contacts["hydrophobic"], COLORS["hydrophobic"], 100)
    _distance_object(prefix + "clashes", contacts["clashes"], COLORS["clash"], 100)
    cmd.zoom(iface_ref, 8)
    return {
        "chain_a_selection": chain_a_sel,
        "chain_b_selection": chain_b_sel,
        "interface_selection": iface_sel,
    }


def contact_rows(kind, contacts):
    rows = []
    for contact in contacts:
        atom_a = contact["atom_a"]
        atom_b = contact["atom_b"]
        rows.append({
            "type": kind,
            "chain_a": atom_a.chain,
            "residue_a": "%s%s" % (atom_a.resn, atom_a.resi),
            "atom_a": atom_a.name,
            "chain_b": atom_b.chain,
            "residue_b": "%s%s" % (atom_b.resn, atom_b.resi),
            "atom_b": atom_b.name,
            "distance": contact["distance"],
            "detail": "%s -- %s" % (_atom_label(atom_a), _atom_label(atom_b)),
        })
    return rows


def add_local_area_to_contact_rows(rows, interface):
    lookup = residue_area_lookup(interface)
    for row in rows:
        area_a = lookup.get((row["chain_a"], row["residue_a"]), 0.0)
        area_b = lookup.get((row["chain_b"], row["residue_b"]), 0.0)
        row["residue_a_dasa"] = area_a
        row["residue_b_dasa"] = area_b
        row["local_pair_area"] = (area_a + area_b) / 2.0
    return rows


def build_pair_table(summary, interface, contacts):
    rows = []
    rows.extend(contact_rows("H-bond", contacts["hbonds"]))
    rows.extend(contact_rows("Salt bridge", contacts["salt_bridges"]))
    rows.extend(contact_rows("Hydrophobic", contacts["hydrophobic"]))
    rows.extend(contact_rows("Close contact", contacts["clashes"]))
    rows = add_local_area_to_contact_rows(rows, interface)
    local_rows = local_area_rows(interface)
    breakdown = summary.get("score_breakdown", {})

    lines = []
    lines.append("[%s] %s chain %s vs chain %s" %
                 (PLUGIN, summary["object"], summary["chain_a"], summary["chain_b"]))
    lines.append("Interface residues: %d" % len(interface["residues"]))
    lines.append("Whole interaction area: %.1f A^2 | Total buried area: %.1f A^2 | Score: %.1f/100" %
                 (summary["interface_area"], summary["total_buried_area"], summary["score"]))
    if breakdown:
        lines.append("Score components: area %.1f + H-bond %.1f + salt %.1f + hydrophobic %.1f - clash %.1f" %
                     (breakdown["area"], breakdown["hbonds"], breakdown["salt_bridges"],
                      breakdown["hydrophobic"], breakdown["clash_penalty"]))
    lines.append("H-bonds: %d | Salt bridges: %d | Hydrophobic contacts: %d | Close contacts: %d" %
                 (len(contacts["hbonds"]), len(contacts["salt_bridges"]),
                  len(contacts["hydrophobic"]), len(contacts["clashes"])))
    lines.append("")
    lines.append("%-14s %-12s %-12s %8s  %s" %
                 ("Type", "Residue A", "Residue B", "Dist(A)", "Atoms"))
    lines.append("-" * 88)
    for row in rows[:250]:
        lines.append("%-14s %-12s %-12s %8.2f  %s" %
                     (row["type"],
                      "%s:%s" % (row["chain_a"], row["residue_a"]),
                      "%s:%s" % (row["chain_b"], row["residue_b"]),
                      row["distance"], row["detail"]))
    if len(rows) > 250:
        lines.append("... %d more interactions omitted; export CSV for all rows." % (len(rows) - 250))

    lines.append("")
    lines.append("Local area by interface residue")
    lines.append("%-8s %-10s %12s %15s %10s" %
                 ("Chain", "Residue", "dASA(A^2)", "LocalArea(A^2)", "Share"))
    lines.append("-" * 62)
    for row in local_rows[:30]:
        lines.append("%-8s %-10s %12.1f %15.1f %9.1f%%" %
                     (row["chain"], row["residue"], row["dasa"],
                      row["local_interface_area"], row["interface_fraction"] * 100.0))

    lines.append("")
    lines.append("Top interface residues by dASA")
    lines.append("%-8s %-10s %10s" % ("Chain", "Residue", "dASA(A^2)"))
    lines.append("-" * 34)
    for residue in sorted(interface["residues"], key=lambda r: r["positive_dasa"], reverse=True)[:30]:
        lines.append("%-8s %-10s %10.1f" %
                     (residue["chain"], "%s%s" % (residue["resn"], residue["resi"]),
                      residue["positive_dasa"]))
    return "\n".join(lines), rows


def export_pair_csv(path, summary, interface, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["summary_key", "summary_value"])
        writer.writerow(["object", summary["object"]])
        writer.writerow(["chain_a", summary["chain_a"]])
        writer.writerow(["chain_b", summary["chain_b"]])
        writer.writerow(["interface_area_A2", "%.2f" % summary["interface_area"]])
        writer.writerow(["total_buried_area_A2", "%.2f" % summary["total_buried_area"]])
        writer.writerow(["score", "%.1f" % summary["score"]])
        if "score_breakdown" in summary:
            writer.writerow(["score_area_component", "%.2f" % summary["score_breakdown"]["area"]])
            writer.writerow(["score_hbond_component", "%.2f" % summary["score_breakdown"]["hbonds"]])
            writer.writerow(["score_salt_component", "%.2f" % summary["score_breakdown"]["salt_bridges"]])
            writer.writerow(["score_hydrophobic_component", "%.2f" % summary["score_breakdown"]["hydrophobic"]])
            writer.writerow(["score_clash_penalty", "%.2f" % summary["score_breakdown"]["clash_penalty"]])
        writer.writerow([])
        writer.writerow(["interaction_type", "chain_a", "residue_a", "atom_a",
                         "chain_b", "residue_b", "atom_b", "distance_A",
                         "residue_a_dASA_A2", "residue_b_dASA_A2",
                         "local_pair_area_A2", "detail"])
        for row in rows:
            writer.writerow([
                row["type"], row["chain_a"], row["residue_a"], row["atom_a"],
                row["chain_b"], row["residue_b"], row["atom_b"],
                "%.2f" % row["distance"],
                "%.2f" % row.get("residue_a_dasa", 0.0),
                "%.2f" % row.get("residue_b_dasa", 0.0),
                "%.2f" % row.get("local_pair_area", 0.0),
                row["detail"],
            ])
        writer.writerow([])
        writer.writerow(["interface_residue_chain", "interface_residue", "dASA_A2",
                         "local_interface_area_A2", "interface_fraction"])
        for residue in local_area_rows(interface):
            writer.writerow([
                residue["chain"], residue["residue"],
                "%.2f" % residue["dasa"],
                "%.2f" % residue["local_interface_area"],
                "%.4f" % residue["interface_fraction"],
            ])


def analyze_pair(complex_obj, chain_a, chain_b, cutoff=1.0, prefix="PPI",
                 visualize=1, export_path="", hbond_cutoff=3.6,
                 salt_cutoff=4.0, hydrophobic_cutoff=4.2, clash_cutoff=2.2,
                 quiet=0):
    pair_prefix = _pair_prefix(prefix, chain_a, chain_b)
    interface = calculate_interface_dasa(complex_obj, chain_a, chain_b, cutoff, pair_prefix)
    contacts = {
        "hbonds": find_hbonds(complex_obj, chain_a, chain_b, hbond_cutoff),
        "salt_bridges": find_salt_bridges(complex_obj, chain_a, chain_b, salt_cutoff),
        "hydrophobic": find_hydrophobic_contacts(complex_obj, chain_a, chain_b, hydrophobic_cutoff),
        "clashes": find_close_contacts(complex_obj, chain_a, chain_b, clash_cutoff),
    }
    summary = {
        "object": complex_obj,
        "chain_a": chain_a,
        "chain_b": chain_b,
        "interface_area": interface["interface_area"],
        "total_buried_area": interface["total_buried_area"],
        "score": score_interface(interface["interface_area"], contacts),
        "score_breakdown": score_breakdown(interface["interface_area"], contacts),
    }
    objects = {}
    if _as_bool(visualize):
        objects = make_visual_objects(complex_obj, chain_a, chain_b, interface, contacts, pair_prefix)
    table, rows = build_pair_table(summary, interface, contacts)
    if export_path:
        export_pair_csv(os.path.expanduser(export_path), summary, interface, rows)
        print("[%s] CSV exported: %s" % (PLUGIN, os.path.expanduser(export_path)))
    if not _as_bool(quiet):
        print(table)
    return {
        "prefix": pair_prefix,
        "summary": summary,
        "interface": interface,
        "contacts": contacts,
        "objects": objects,
        "rows": rows,
        "table": table,
    }


def build_summary_table(complex_obj, chains, rows):
    total_area = sum(row["interface_area"] for row in rows)
    total_buried = sum(row["buried_area"] for row in rows)
    lines = []
    lines.append("[%s] Pairwise interaction summary for %s chains %s" %
                 (PLUGIN, complex_obj, "+".join(chains)))
    lines.append("Total pairwise interaction area: %.1f A^2 | Total buried area: %.1f A^2" %
                 (total_area, total_buried))
    lines.append("%-8s %8s %10s %7s %6s %12s %8s %7s %7s %7s  %s" %
                 ("Pair", "Residues", "Area(A2)", "HBond", "Salt",
                  "Hydrophobic", "Clash", "Score", "AreaSc", "PolarSc", "Objects"))
    lines.append("-" * 118)
    for row in rows:
        lines.append("%-8s %8d %10.1f %7d %6d %12d %8d %7.1f %7.1f %7.1f  %s*" %
                     (row["pair"], row["interface_residues"], row["interface_area"],
                      row["hbonds"], row["salt_bridges"], row["hydrophobic"],
                      row["clashes"], row["score"], row.get("score_area", 0.0),
                      row.get("score_polar", 0.0), row["prefix"]))
    lines.extend(build_local_hotspot_lines(rows, include_object=False))
    return "\n".join(lines)


def _format_hotspot_rows(hotspots, limit=5):
    if not hotspots:
        return "none"
    parts = []
    for item in hotspots[:limit]:
        parts.append("%s:%s %.1fA2" %
                     (item["chain"], item["residue"], item["local_interface_area"]))
    return "; ".join(parts)


def build_local_hotspot_lines(rows, include_object=True):
    lines = []
    lines.append("")
    lines.append("Local interface area hot spots (residue-level LocalArea = dASA/2)")
    if include_object:
        lines.append("%-12s %-8s  %s" % ("Object", "Pair", "Top local residues"))
        lines.append("-" * 90)
        for row in rows:
            lines.append("%-12s %-8s  %s" %
                         (row.get("object", ""), row["pair"],
                          _format_hotspot_rows(row.get("local_hotspots", []))))
    else:
        lines.append("%-8s  %s" % ("Pair", "Top local residues"))
        lines.append("-" * 78)
        for row in rows:
            lines.append("%-8s  %s" %
                         (row["pair"], _format_hotspot_rows(row.get("local_hotspots", []))))
    return lines


def build_batch_table(rows):
    total_area = sum(row["interface_area"] for row in rows)
    total_buried = sum(row["buried_area"] for row in rows)
    lines = []
    lines.append("[%s] Batch interaction summary" % PLUGIN)
    lines.append("Total interaction area across all jobs: %.1f A^2 | Total buried area: %.1f A^2" %
                 (total_area, total_buried))
    lines.append("%-12s %-8s %8s %10s %7s %6s %12s %8s %7s %7s %7s  %s" %
                 ("Object", "Pair", "Residues", "Area(A2)", "HBond", "Salt",
                  "Hydrophobic", "Clash", "Score", "AreaSc", "PolarSc", "Objects"))
    lines.append("-" * 132)
    for row in rows:
        lines.append("%-12s %-8s %8d %10.1f %7d %6d %12d %8d %7.1f %7.1f %7.1f  %s*" %
                     (row["object"], row["pair"], row["interface_residues"],
                      row["interface_area"], row["hbonds"], row["salt_bridges"],
                      row["hydrophobic"], row["clashes"], row["score"],
                      row.get("score_area", 0.0), row.get("score_polar", 0.0),
                      row["prefix"]))
    lines.extend(build_local_hotspot_lines(rows, include_object=True))
    return "\n".join(lines)


def export_summary_csv(path, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["pair", "interface_residues", "interface_area_A2",
                         "buried_area_A2", "hbonds", "salt_bridges",
                         "hydrophobic_contacts", "close_contacts", "score",
                         "score_area", "score_polar", "score_hydrophobic",
                         "score_clash_penalty",
                         "top_local_area_residues",
                         "object_prefix"])
        for row in rows:
            writer.writerow([
                row["pair"], row["interface_residues"],
                "%.2f" % row["interface_area"], "%.2f" % row["buried_area"],
                row["hbonds"], row["salt_bridges"], row["hydrophobic"],
                row["clashes"], "%.1f" % row["score"],
                "%.2f" % row.get("score_area", 0.0),
                "%.2f" % row.get("score_polar", 0.0),
                "%.2f" % row.get("score_hydrophobic", 0.0),
                "%.2f" % row.get("score_clash_penalty", 0.0),
                _format_hotspot_rows(row.get("local_hotspots", [])),
                row["prefix"],
            ])


def export_batch_summary_csv(path, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["object", "pair", "interface_residues", "interface_area_A2",
                         "buried_area_A2", "hbonds", "salt_bridges",
                         "hydrophobic_contacts", "close_contacts", "score",
                         "score_area", "score_polar", "score_hydrophobic",
                         "score_clash_penalty",
                         "top_local_area_residues",
                         "object_prefix"])
        for row in rows:
            writer.writerow([
                row["object"], row["pair"], row["interface_residues"],
                "%.2f" % row["interface_area"], "%.2f" % row["buried_area"],
                row["hbonds"], row["salt_bridges"], row["hydrophobic"],
                row["clashes"], "%.1f" % row["score"],
                "%.2f" % row.get("score_area", 0.0),
                "%.2f" % row.get("score_polar", 0.0),
                "%.2f" % row.get("score_hydrophobic", 0.0),
                "%.2f" % row.get("score_clash_penalty", 0.0),
                _format_hotspot_rows(row.get("local_hotspots", [])),
                row["prefix"],
            ])


def ppi_analyze(complex_obj="all", chains="A+B", cutoff=1.0, prefix="PPI",
                visualize=1, export_path="", export_dir="", summary_csv="",
                hbond_cutoff=3.6, salt_cutoff=4.0,
                hydrophobic_cutoff=4.2, clash_cutoff=2.2):
    """Analyze one chain pair or all pairwise interactions among several chains."""
    if isinstance(complex_obj, str) and "," in complex_obj:
        complex_obj, chains, inline_options = _parse_inline_command(complex_obj, chains)
        cutoff = _float_option(inline_options, "cutoff", cutoff)
        prefix = _string_option(inline_options, "prefix", prefix)
        visualize = _string_option(inline_options, "visualize", visualize)
        export_path = _string_option(inline_options, "export_path", export_path)
        export_dir = _string_option(inline_options, "export_dir", export_dir)
        summary_csv = _string_option(inline_options, "summary_csv", summary_csv)
        hbond_cutoff = _float_option(inline_options, "hbond_cutoff", hbond_cutoff)
        salt_cutoff = _float_option(inline_options, "salt_cutoff", salt_cutoff)
        hydrophobic_cutoff = _float_option(inline_options, "hydrophobic_cutoff", hydrophobic_cutoff)
        clash_cutoff = _float_option(inline_options, "clash_cutoff", clash_cutoff)

    chain_values = _parse_chains(chains)
    if len(chain_values) == 2:
        return analyze_pair(
            complex_obj, chain_values[0], chain_values[1],
            cutoff=cutoff, prefix=prefix, visualize=visualize,
            export_path=export_path, hbond_cutoff=hbond_cutoff,
            salt_cutoff=salt_cutoff, hydrophobic_cutoff=hydrophobic_cutoff,
            clash_cutoff=clash_cutoff,
        )

    if export_dir:
        export_dir = os.path.expanduser(export_dir)
        if not os.path.isdir(export_dir):
            os.makedirs(export_dir)

    results = {}
    summary_rows = []
    for chain_a, chain_b in itertools.combinations(chain_values, 2):
        pair_prefix = _pair_prefix(prefix, chain_a, chain_b)
        pair_export = os.path.join(export_dir, pair_prefix + ".csv") if export_dir else ""
        result = analyze_pair(
            complex_obj, chain_a, chain_b,
            cutoff=cutoff, prefix=prefix, visualize=visualize,
            export_path=pair_export, hbond_cutoff=hbond_cutoff,
            salt_cutoff=salt_cutoff, hydrophobic_cutoff=hydrophobic_cutoff,
            clash_cutoff=clash_cutoff, quiet=1,
        )
        row = {
            "pair": "%s-%s" % (chain_a, chain_b),
            "prefix": pair_prefix,
            "interface_residues": len(result["interface"]["residues"]),
            "interface_area": result["summary"]["interface_area"],
            "buried_area": result["summary"]["total_buried_area"],
            "hbonds": len(result["contacts"]["hbonds"]),
            "salt_bridges": len(result["contacts"]["salt_bridges"]),
            "hydrophobic": len(result["contacts"]["hydrophobic"]),
            "clashes": len(result["contacts"]["clashes"]),
            "score": result["summary"]["score"],
            "score_area": result["summary"]["score_breakdown"]["area"],
            "score_polar": (
                result["summary"]["score_breakdown"]["hbonds"] +
                result["summary"]["score_breakdown"]["salt_bridges"]
            ),
            "score_hydrophobic": result["summary"]["score_breakdown"]["hydrophobic"],
            "score_clash_penalty": result["summary"]["score_breakdown"]["clash_penalty"],
            "local_hotspots": local_area_rows(result["interface"])[:5],
        }
        summary_rows.append(row)
        results[row["pair"]] = result

    summary_rows.sort(key=lambda r: r["score"], reverse=True)
    table = build_summary_table(complex_obj, chain_values, summary_rows)
    print(table)
    if summary_csv:
        export_summary_csv(os.path.expanduser(summary_csv), summary_rows)
        print("[%s] Summary CSV exported: %s" % (PLUGIN, os.path.expanduser(summary_csv)))
    return {
        "chains": chain_values,
        "summary_rows": summary_rows,
        "results": results,
        "table": table,
    }


def ppi_batch(specs="", cutoff=1.0, prefix="PPI", visualize=1,
              export_dir="", summary_csv="", hbond_cutoff=3.6,
              salt_cutoff=4.0, hydrophobic_cutoff=4.2, clash_cutoff=2.2):
    """Analyze two selected chains across multiple loaded objects/PDBs.

    Example:
        ppi_batch "001:A+B;002:A+B;003:C+D"
    """
    jobs = _parse_batch_specs(specs)
    if export_dir:
        export_dir = os.path.expanduser(export_dir)
        if not os.path.isdir(export_dir):
            os.makedirs(export_dir)

    rows = []
    results = {}
    for obj, chains in jobs:
        chain_values = _parse_chains(chains)
        if len(chain_values) != 2:
            raise ValueError("Batch item must contain exactly two chains: %s:%s" % (obj, chains))
        chain_a, chain_b = chain_values
        job_prefix = _safe_name("%s%s%s%s" % (prefix, obj, chain_a, chain_b), "PPI")
        pair_export = os.path.join(export_dir, job_prefix + ".csv") if export_dir else ""
        result = analyze_pair(
            obj, chain_a, chain_b, cutoff=cutoff, prefix=job_prefix,
            visualize=visualize, export_path=pair_export,
            hbond_cutoff=hbond_cutoff, salt_cutoff=salt_cutoff,
            hydrophobic_cutoff=hydrophobic_cutoff,
            clash_cutoff=clash_cutoff, quiet=1,
        )
        row = {
            "object": obj,
            "pair": "%s-%s" % (chain_a, chain_b),
            "prefix": result["prefix"],
            "interface_residues": len(result["interface"]["residues"]),
            "interface_area": result["summary"]["interface_area"],
            "buried_area": result["summary"]["total_buried_area"],
            "hbonds": len(result["contacts"]["hbonds"]),
            "salt_bridges": len(result["contacts"]["salt_bridges"]),
            "hydrophobic": len(result["contacts"]["hydrophobic"]),
            "clashes": len(result["contacts"]["clashes"]),
            "score": result["summary"]["score"],
            "score_area": result["summary"]["score_breakdown"]["area"],
            "score_polar": (
                result["summary"]["score_breakdown"]["hbonds"] +
                result["summary"]["score_breakdown"]["salt_bridges"]
            ),
            "score_hydrophobic": result["summary"]["score_breakdown"]["hydrophobic"],
            "score_clash_penalty": result["summary"]["score_breakdown"]["clash_penalty"],
            "local_hotspots": local_area_rows(result["interface"])[:5],
        }
        rows.append(row)
        results["%s:%s" % (obj, row["pair"])] = result

    rows.sort(key=lambda r: (r["object"], r["pair"]))
    table = build_batch_table(rows)
    print(table)
    if summary_csv:
        export_batch_summary_csv(os.path.expanduser(summary_csv), rows)
        print("[%s] Batch summary CSV exported: %s" % (PLUGIN, os.path.expanduser(summary_csv)))
    return {"jobs": jobs, "summary_rows": rows, "results": results, "table": table}


def ppi_score_help():
    text = """
[Protein Interaction Analyzer] How to compare interface strength

Use Score for quick ranking, then inspect the components:
  Score        Overall heuristic score, 0-100.
  AreaSc       Contribution from whole interaction area, max 40.
  PolarSc      H-bond + salt bridge contribution, max 36.
  Hydrophobic  Hydrophobic-contact contribution, max 20.
  ClashPenalty Penalty from very close contacts, max 30.

Practical interpretation:
  Stronger-looking interface: high Score, high area, many H-bonds/salt bridges,
  enough hydrophobic contacts, and low clash penalty.

Local area:
  dASA(A^2) is residue-level buried surface area.
  LocalArea(A^2) is dASA/2, used as a residue-level estimate of contribution
  to the interface area. Large LocalArea residues are likely interface hot spots
  worth inspecting or mutating.

Important:
  This is a structure-based heuristic, not a binding free energy calculation.
  Compare related models/complexes under the same cutoffs.
"""
    print(text)
    return text


def _import_tk():
    try:
        import tkinter as tk
        from tkinter import ttk, filedialog, messagebox
        return tk, ttk, filedialog, messagebox, None
    except Exception as exc:
        first_error = exc
    try:
        import Tkinter as tk
        import ttk
        import tkFileDialog as filedialog
        import tkMessageBox as messagebox
        return tk, ttk, filedialog, messagebox, None
    except Exception:
        return None, None, None, None, first_error


class PPIAnalyzerPanel(object):
    def __init__(self, parent=None):
        tk, ttk, filedialog, messagebox, error = _import_tk()
        if error is not None:
            raise RuntimeError("Tk/Tcl is unavailable in this PyMOL build: %s" % error)
        self.tk = tk
        self.filedialog = filedialog
        self.messagebox = messagebox
        self.last_result = None

        self.top = tk.Toplevel(parent)
        self.top.title(PLUGIN)
        self.top.geometry("820x720")
        self.top.columnconfigure(1, weight=1)
        row = 0

        ttk.Label(self.top, text="PDB ID / local file").grid(row=row, column=0, sticky="w", padx=10, pady=5)
        load_frame = ttk.Frame(self.top)
        load_frame.grid(row=row, column=1, sticky="ew", padx=10, pady=5)
        load_frame.columnconfigure(0, weight=1)
        self.load_var = tk.StringVar()
        ttk.Entry(load_frame, textvariable=self.load_var).grid(row=0, column=0, sticky="ew")
        ttk.Button(load_frame, text="Fetch/Load", command=self._load_structure).grid(row=0, column=1, padx=(6, 0))
        ttk.Button(load_frame, text="Browse", command=self._browse).grid(row=0, column=2, padx=(6, 0))
        row += 1

        ttk.Label(self.top, text="Object").grid(row=row, column=0, sticky="w", padx=10, pady=5)
        self.object_var = tk.StringVar()
        self.object_box = ttk.Combobox(self.top, textvariable=self.object_var, values=self._objects())
        self.object_box.grid(row=row, column=1, sticky="ew", padx=10, pady=5)
        objects = self._objects()
        if objects:
            self.object_var.set(objects[0])
        self.object_box.bind("<<ComboboxSelected>>", self._refresh_chains)
        row += 1

        ttk.Label(self.top, text="Chains").grid(row=row, column=0, sticky="w", padx=10, pady=5)
        chain_frame = ttk.Frame(self.top)
        chain_frame.grid(row=row, column=1, sticky="ew", padx=10, pady=5)
        chain_frame.columnconfigure(0, weight=1)
        self.chains_var = tk.StringVar(value=self._default_chains())
        ttk.Entry(chain_frame, textvariable=self.chains_var).grid(row=0, column=0, sticky="ew")
        ttk.Button(chain_frame, text="Refresh", command=self._refresh_chains).grid(row=0, column=1, padx=(6, 0))
        row += 1

        self.vars = {}
        fields = [
            ("dASA cutoff (A^2)", "cutoff", "1.0"),
            ("H-bond cutoff (A)", "hbond", "3.6"),
            ("Salt bridge cutoff (A)", "salt", "4.0"),
            ("Hydrophobic cutoff (A)", "hydrophobic", "4.2"),
            ("Close contact cutoff (A)", "clash", "2.2"),
            ("Output prefix", "prefix", "PPI"),
        ]
        for label, key, default in fields:
            ttk.Label(self.top, text=label).grid(row=row, column=0, sticky="w", padx=10, pady=5)
            self.vars[key] = tk.StringVar(value=default)
            ttk.Entry(self.top, textvariable=self.vars[key]).grid(row=row, column=1, sticky="ew", padx=10, pady=5)
            row += 1

        self.visualize_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            self.top,
            text="Create chain/interface selections, contact objects, and labels in PyMOL",
            variable=self.visualize_var,
        ).grid(row=row, column=0, columnspan=2, sticky="w", padx=10, pady=5)
        row += 1

        buttons = ttk.Frame(self.top)
        buttons.grid(row=row, column=0, columnspan=2, sticky="ew", padx=10, pady=8)
        buttons.columnconfigure(0, weight=1)
        buttons.columnconfigure(1, weight=1)
        ttk.Button(buttons, text="Analyze", command=self._run).grid(row=0, column=0, sticky="ew", padx=(0, 5))
        ttk.Button(buttons, text="Export table", command=self._export).grid(row=0, column=1, sticky="ew", padx=(5, 0))
        row += 1

        self.output = tk.Text(self.top, height=26, wrap="none")
        self.output.grid(row=row, column=0, columnspan=2, sticky="nsew", padx=10, pady=(0, 10))
        self.top.rowconfigure(row, weight=1)

    def _objects(self):
        return cmd.get_names("objects") or []

    def _chains(self):
        obj = self.object_var.get()
        if not obj and self._objects():
            obj = self._objects()[0]
        try:
            return [c for c in cmd.get_chains(obj) if c]
        except Exception:
            return []

    def _default_chains(self):
        chains = self._chains()
        if len(chains) >= 2:
            return "+".join(chains[:2])
        return "A+B"

    def _refresh_chains(self, _event=None):
        objects = self._objects()
        self.object_box.configure(values=objects)
        if not self.object_var.get() and objects:
            self.object_var.set(objects[0])
        self.chains_var.set(self._default_chains())

    def _browse(self):
        path = self.filedialog.askopenfilename(
            title="Open structure",
            filetypes=[("Structure files", "*.pdb *.cif *.mmcif"), ("All files", "*.*")]
        )
        if path:
            self.load_var.set(path)

    def _load_structure(self):
        target = self.load_var.get().strip()
        if not target:
            return
        try:
            expanded = os.path.expanduser(target)
            if os.path.exists(expanded):
                obj = ppi_load(expanded)
            else:
                obj = ppi_fetch(target)
            self.object_var.set(obj)
            self._refresh_chains()
        except Exception as exc:
            self.messagebox.showerror(PLUGIN, str(exc))

    def _options(self):
        return {
            "cutoff": float(self.vars["cutoff"].get()),
            "prefix": self.vars["prefix"].get(),
            "visualize": self.visualize_var.get(),
            "hbond_cutoff": float(self.vars["hbond"].get()),
            "salt_cutoff": float(self.vars["salt"].get()),
            "hydrophobic_cutoff": float(self.vars["hydrophobic"].get()),
            "clash_cutoff": float(self.vars["clash"].get()),
        }

    def _run(self):
        try:
            result = ppi_analyze(self.object_var.get(), self.chains_var.get(), **self._options())
            self.last_result = result
            self.output.delete("1.0", self.tk.END)
            self.output.insert(self.tk.END, result["table"])
        except Exception as exc:
            self.messagebox.showerror(PLUGIN, str(exc))

    def _export(self):
        if not self.last_result:
            return
        path = self.filedialog.asksaveasfilename(
            title="Export table",
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv"), ("All files", "*.*")]
        )
        if not path:
            return
        try:
            if "summary_rows" in self.last_result:
                export_summary_csv(path, self.last_result["summary_rows"])
            else:
                export_pair_csv(
                    path,
                    self.last_result["summary"],
                    self.last_result["interface"],
                    self.last_result["rows"],
                )
            self.messagebox.showinfo(PLUGIN, "Exported: %s" % path)
        except Exception as exc:
            self.messagebox.showerror(PLUGIN, str(exc))


def ppi_panel():
    _tk, _ttk, _filedialog, _messagebox, error = _import_tk()
    if error is not None:
        print("[%s] GUI panel cannot open because Tk/Tcl is unavailable: %s" % (PLUGIN, error))
        print("Use commands instead:")
        print("  ppi_fetch 1brs")
        print("  ppi_analyze 1brs, A+B")
        print("  ppi_analyze 1brs, A+B+C+D")
        print('  ppi_batch "001:A+B;002:A+B;003:C+D"')
        return None
    parent = None
    try:
        from pymol import plugins
        parent = plugins.get_tk_root()
    except Exception:
        pass
    return PPIAnalyzerPanel(parent)


def __init_plugin__(app=None):
    cmd.extend("ppi_panel", ppi_panel)
    cmd.extend("ppi_fetch", ppi_fetch)
    cmd.extend("ppi_load", ppi_load)
    cmd.extend("ppi_analyze", ppi_analyze)
    cmd.extend("ppi_batch", ppi_batch)
    cmd.extend("ppi_score_help", ppi_score_help)
    try:
        app.menuBar.addmenuitem("Plugin", "command", PLUGIN, label=PLUGIN, command=ppi_panel)
    except Exception:
        pass


cmd.extend("ppi_panel", ppi_panel)
cmd.extend("ppi_fetch", ppi_fetch)
cmd.extend("ppi_load", ppi_load)
cmd.extend("ppi_analyze", ppi_analyze)
cmd.extend("ppi_batch", ppi_batch)
cmd.extend("ppi_score_help", ppi_score_help)
