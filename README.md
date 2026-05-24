# Protein Interaction Analyzer

Copyright (c) yangyu. All rights reserved.

Protein Interaction Analyzer is a PyMOL plugin for quick inspection and comparison of protein-protein interfaces.

## Features

- Interface residue detection based on dASA
- Whole interaction area and total buried area
- Residue-level local interface area hot spots
- Hydrogen bond detection
- Salt bridge / ionic interaction detection
- Hydrophobic contact detection
- Close contact / potential clash detection
- Heuristic 0-100 interface strength score
- Pairwise analysis among multiple chains in one structure
- Batch comparison across multiple loaded PyMOL objects
- Region-based interaction search from any user-defined selection
- PTM-site interaction search for common modified residues such as SEP/TPO/PTR
- CSV export

## Recommended Plugin File

Use the latest compatibility version:

```text
ProteinInteractionPlugin_v3.py
```

Older versions are kept for traceability:

- `ProteinInteractionPlugin.py`: original PyMOL version
- `ProteinInteractionPlugin_v2.py`: compatibility fix for PyMOL Atom objects without `.model`
- `ProteinInteractionPlugin_v3.py`: v2 plus more robust interface residue labels

## Quick Start

Load the plugin in PyMOL:

```pymol
run /path/to/interfacefinder/ProteinInteractionPlugin_v3.py
```

Fetch a test structure:

```pymol
ppi_fetch 1brs
```

Analyze one interface:

```pymol
ppi_analyze 1brs, A+D
```

Batch analysis:

```pymol
ppi_batch "001:A+B;002:A+B;003:C+D"
```

Search interactions around a selected region:

```pymol
select my_region, 1brs and chain A and resi 35-45
ppi_region my_region
```

Search interactions around PTM sites:

```pymol
ppi_ptm 7mp9
ppi_ptm "object=7mp9; ptm=phospho"
```

## Documentation

Read the full user guide:

```text
USER_GUIDE.md
```

Version notes:

```text
ProteinInteractionPlugin_v2_CompatibilityNote.md
ProteinInteractionPlugin_v3_ReleaseNote.md
```

## Important Notes

The score is a geometry-based heuristic for comparing related structures or models. It is not a binding free-energy calculation and should not be interpreted as experimental affinity.
