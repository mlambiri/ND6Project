# Simulation helpers (membrane sizing / patch build)

These helpers are for building a **lipid bilayer patch** sized to a Complex I system for MD/NAMD-style simulations.

## Requirement implemented

The membrane patch dimensions are computed in the **inferred membrane plane**:

- infer the membrane normal from **PCA of CA atoms in ND1..ND6**
- compute ND extents in that plane
- `membrane_x = (extent_u + 50 A)`
- `membrane_y = (extent_v + 50 A)`

This corresponds to a **+25 A margin on each side** (default).

## 1) Compute required patch size (no VMD needed)

```bash
python simulation/calc_membrane_patch_size.py --pdb /path/to/system.pdb --margin 25
```

This prints the ND-plane extents and the recommended membrane patch X/Y.

Recommended inputs (this repo's **full Complex I** models):

- WT (protein-only): `output/playwright/chatgpt_botprompts/models/complexI_9TI4_WT_heavy_proteinOnly.pdb`
- WT (full heavy; includes native lipids/cofactors): `output/playwright/chatgpt_botprompts/models/complexI_9TI4_WT_heavy.pdb`
- Mutants (full heavy): `output/playwright/chatgpt_botprompts/models/complexI_9TI4_ND1_A52T_heavy.pdb`, `output/playwright/chatgpt_botprompts/models/complexI_9TI4_ND4_R340H_heavy.pdb`, `output/playwright/chatgpt_botprompts/models/complexI_9TI4_ND6_M64V_heavy.pdb`

Example (full Complex I, WT):

```bash
python simulation/calc_membrane_patch_size.py --pdb output/playwright/chatgpt_botprompts/models/complexI_9TI4_WT_heavy_proteinOnly.pdb --margin 25
```

## 2) Build a membrane patch in VMD (membrane plugin)

Requires VMD with the `membrane` plugin available.

1) Run step (1) and note the recommended `x=` and `y=` values.

2) Build the patch using explicit sizes (recommended):

```bash
vmd -dispdev text -e simulation/vmd_build_membrane_patch.tcl -args --pdb output/playwright/chatgpt_botprompts/models/complexI_9TI4_WT_heavy_proteinOnly.pdb --x <X_A> --y <Y_A> --lipid POPC --out simulation/out/complexI_9TI4_membrane --place 0
```

Outputs (prefix = `--out`):

- `<out>.psf`
- `<out>.pdb`
- `<out>_placed.pdb` (optionally translated to be centered on the protein in X/Y)

Notes:

- Assumes the **membrane plane is XY** and the **normal is Z** in the input coordinates (typical NAMD setups).
- The default extent selection excludes common **water/ion/lipid** residue names; pass `--sel` to override.
- This script **only builds the bilayer patch**; embedding + resolvation/ionization is intentionally left to your system-build workflow (CHARMM-GUI, VMD solvate/autoionize, etc.).

## 3) Place the patch into the Complex I coordinate frame (Python)

The membrane patch produced by the VMD `membrane` plugin is centered near `z~0` and aligned to XYZ axes. Your Complex I models are not, so you generally need to **rotate + translate** the patch to match the protein coordinate frame.

This script:

- infers the membrane plane from **ND CA atoms**
- rotates the patch so its normal matches that plane
- centers the patch on the **ND bounding-box center** in the membrane plane (so margins apply on all sides)

Example:

```bash
python simulation/place_membrane_patch.py --protein-pdb output/playwright/chatgpt_botprompts/models/complexI_9TI4_WT_heavy_proteinOnly.pdb --patch-pdb simulation/out/complexI_9TI4_membrane.pdb --out-pdb simulation/out/complexI_9TI4_membrane_placed.pdb
```

## 4) (Optional) Write a lipids-only patch (strip waters) in VMD

If you want a smaller membrane file for visualization (no waters), you can export a selection to a new PSF+PDB:

```bash
vmd -dispdev text -e simulation/vmd_export_selection_psf_pdb.tcl -args --psf simulation/out/complexI_9TI4_membrane.psf --pdb simulation/out/complexI_9TI4_membrane_placed.pdb --sel "not water" --out simulation/out/complexI_9TI4_membrane_placed_lipidsOnly
```
