# Simulation helpers (membrane sizing / patch build)

These helpers are for building a **lipid bilayer patch** sized to a Complex I system for MD/NAMD-style simulations.

## Requirement implemented

The membrane patch X/Y dimensions are computed from the Complex I bounding box in the membrane plane:

- `membrane_x = (complex_x + 50 A)`
- `membrane_y = (complex_y + 50 A)`

This corresponds to a **+25 A margin on each side** (default).

## 1) Compute required patch size (no VMD needed)

```bash
python simulation/calc_membrane_patch_size.py --pdb /path/to/system.pdb --margin 25
```

This prints the Complex I (non-water / non-ion / non-lipid) extents and the recommended membrane X/Y.

Recommended inputs (this repo’s **full Complex I** models):

- WT (protein-only): `output/playwright/chatgpt_botprompts/models/complexI_9TI4_WT_heavy_proteinOnly.pdb`
- WT (full heavy; includes native lipids/cofactors): `output/playwright/chatgpt_botprompts/models/complexI_9TI4_WT_heavy.pdb`
- Mutants (full heavy): `output/playwright/chatgpt_botprompts/models/complexI_9TI4_ND1_A52T_heavy.pdb`, `output/playwright/chatgpt_botprompts/models/complexI_9TI4_ND4_R340H_heavy.pdb`, `output/playwright/chatgpt_botprompts/models/complexI_9TI4_ND6_M64V_heavy.pdb`

Example (full Complex I, WT):

```bash
python simulation/calc_membrane_patch_size.py --pdb output/playwright/chatgpt_botprompts/models/complexI_9TI4_WT_heavy_proteinOnly.pdb --margin 25
```

## 2) Build a membrane patch in VMD (membrane plugin)

Requires VMD with the `membrane` plugin available.

```bash
vmd -dispdev text -e simulation/vmd_build_membrane_patch.tcl -args --pdb output/playwright/chatgpt_botprompts/models/complexI_9TI4_WT_heavy_proteinOnly.pdb --margin 25 --lipid POPC --out simulation/out/complexI_9TI4_membrane
```

Outputs (prefix = `--out`):

- `<out>.psf`
- `<out>.pdb`
- `<out>_placed.pdb` (optionally translated to be centered on the protein in X/Y)

Notes:

- Assumes the **membrane plane is XY** and the **normal is Z** in the input coordinates (typical NAMD setups).
- The default extent selection excludes common **water/ion/lipid** residue names; pass `--sel` to override.
- This script **only builds the bilayer patch**; embedding + resolvation/ionization is intentionally left to your system-build workflow (CHARMM-GUI, VMD solvate/autoionize, etc.).
