# Complex I (human) + ND6 M64V starter models

Exported from [RCSB PDB 9TI4](https://www.rcsb.org/structure/9TI4)

This folder contains "clean" (heavy-atoms-only) coordinate files exported from RCSB **PDB 9TI4** (human mitochondrial Complex I; cryo-EM) and a simple ND6 mutation helper.

## Files

- `complexI_9TI4_WT_heavy.pdb`: full Complex I, heavy atoms only (HOH removed)
- `complexI_9TI4_ND1_A52T_heavy.pdb`: same, but **MT-ND1 ALA52->THR** applied
- `complexI_9TI4_ND6_M64V_heavy.pdb`: same, but **MT-ND6 MET64→VAL** applied
- `complexI_9TI4_WT_heavy_proteinOnly.pdb`: protein-only (ATOM records), heavy atoms only
- `complexI_9TI4_ND1_A52T_heavy_proteinOnly.pdb`: protein-only (ATOM records), heavy atoms only, mutated
- `complexI_9TI4_ND6_M64V_heavy_proteinOnly.pdb`: protein-only (ATOM records), heavy atoms only, mutated
- `nd1_chain_s_WT_heavy.pdb`: ND1 only (chain `s`)
- `nd1_chain_s_A52T_heavy.pdb`: ND1 only (chain `s`), mutated
- `nd6_chain_m_WT_heavy.pdb`: ND6 only (chain `m`)
- `nd6_chain_m_M64V_heavy.pdb`: ND6 only (chain `m`), mutated
- `complexI_9TI4_chain_map.json`: author chain IDs → PDB chain IDs (identity for 9TI4)

## Variant mapping

LHON mtDNA variant **m.14484T>C** in **MT-ND6** maps to protein change **p.Met64Val (M64V)**.

In **9TI4**, ND6 corresponds to **chain `m`**, and the mutated residue is **resid 64**.

LHON mtDNA variant **m.3460G>A** in **MT-ND1** maps to protein change **p.Ala52Thr (A52T)**.

In **9TI4**, ND1 corresponds to **chain `s`**, and the mutated residue is **resid 52**.

## Regenerating these outputs

From the repo root:

```bash
python3 output/playwright/chatgpt_botprompts/models/build_complexI_9TI4_models.py
```

## VMD/CHARMM rebuild (optional)

The PDBs here have no hydrogens; you generally want to rebuild/protonate with your forcefield.

For ND6-only, run (from this folder) in VMD:

```bash
export CHARMM_TOP=/path/to/top_all36_prot.rtf
vmd -dispdev text -e vmd_psfgen_mutate_nd6_M64V.tcl
```

This will write `nd6_chain_m_M64V.psf` and `nd6_chain_m_M64V_rebuilt.pdb`.

For ND1-only, run (from this folder) in VMD:

```bash
export CHARMM_TOP=/path/to/top_all36_prot.rtf
vmd -dispdev text -e vmd_psfgen_mutate_nd1_A52T.tcl
```

This will write `nd1_chain_s_A52T.psf` and `nd1_chain_s_A52T_rebuilt.pdb`.

## VMD visualization

From this folder:

```bash
vmd -e vmd_view_complexI_9TI4_WT.tcl
vmd -e vmd_view_complexI_9TI4_ND6_M64V.tcl
vmd -e vmd_compare_complexI_9TI4_WT_vs_M64V.tcl
vmd -e vmd_view_complexI_9TI4_ND1_A52T.tcl
vmd -e vmd_compare_complexI_9TI4_WT_vs_ND1_A52T.tcl
vmd -e vmd_view_complexI_9TI4_features.tcl
vmd -e vmd_compare_complexI_9TI4_WT_vs_M64V_features.tcl
```

By default, the basic scripts (`vmd_view_*` / `vmd_compare_*`) load the `*_proteinOnly.pdb` files (less noisy in VMD). The `*_features.tcl` scripts default to the full heavy-atom models so lipids/cofactors are visible. Pass explicit PDB paths via `-args` as needed.

## VMD: `vmd_view_complexI_9TI4_basic_003.tcl` (ND1–ND6 + membrane toggles)

This script is a label-free “basic” view that lets you turn on/off:

- **ND subunits** (`ND1`…`ND6`)
- **Membrane/peripheral arm** approximation (based on proximity to lipids, or an ND-neighborhood fallback for protein-only models)
- **Membrane components** (native lipids) and **cofactors / Fe-S clusters**

### Run (Windows)

From this folder (`...\\output\\playwright\\chatgpt_botprompts\\models`):

```powershell
cd C:\home\workspace\8196399\output\playwright\chatgpt_botprompts\models
& "C:\Program Files\University of Illinois\VMD2\vmd.exe" -dispdev win -e vmd_view_complexI_9TI4_basic_003.tcl
```

### Command-line options

All options go **after** `-args`.

Common examples:

```powershell
# ND6 only + lipids (good for “embedded in membrane” view)
& "C:\Program Files\University of Illinois\VMD2\vmd.exe" -dispdev win -e vmd_view_complexI_9TI4_basic_003.tcl -args --nd ND6 --show nd,lipids

# ND4 only + membrane arm + lipids (hide cofactors/FeS clutter)
& "C:\Program Files\University of Illinois\VMD2\vmd.exe" -dispdev win -e vmd_view_complexI_9TI4_basic_003.tcl -args --nd ND4 --arm membrane --show nd,arm,lipids --hide cofactors,fes

# All ND subunits + both arms + lipids
& "C:\Program Files\University of Illinois\VMD2\vmd.exe" -dispdev win -e vmd_view_complexI_9TI4_basic_003.tcl -args --nd all --show all --arm both

# Explicitly choose which PDB to load
& "C:\Program Files\University of Illinois\VMD2\vmd.exe" -dispdev win -e vmd_view_complexI_9TI4_basic_003.tcl -args complexI_9TI4_WT_heavy.pdb --show all --arm both
```

Flags summary:

- `--nd ND6` (or `ND1`…`ND5`, `all`, `none`)
- `--show default|all|context,nd,arm,lipids,cofactors,fes`
- `--hide lipids` (comma-separated list)
- `--arm membrane|peripheral|both|off`
- `--arm-include-nd 1` (include ND chains inside the arm representation; default `0`)
- `--membrane-dist 5.0` (distance in Å for “membrane arm” proximity selection)

Note: lipids/cofactors are only present in the full `*_heavy.pdb` models (not in `*_proteinOnly.pdb`).

### Interactive toggles (Tk Console)

After the script loads, open **Extensions → Tk Console** and use:

```tcl
ci_list

# ND subunits
ci_focus_nd ND6      ;# show only ND6
ci_show ND4
ci_show nd           ;# show all ND1..ND6
ci_hide nd

# Membrane / arm
ci_show lipids
ci_hide lipids
ci_show arm_membrane
ci_show arm_peripheral
ci_hide arm

# Cofactors / Fe-S
ci_show cofactors
ci_hide cofactors
ci_show fes
ci_hide fes
```

### Same toggles in other scripts

Most visualization scripts in this folder now share the same CLI flags (`--show/--hide/--nd/--arm`) and Tk Console commands (`ci_show`, `ci_hide`, `ci_focus_nd`) via `ci_showhide.tcl`.

Compatibility wrappers:

- `vmd_view_complexI_9TI4_basic.tcl` and `vmd_view_complexI_9TI4_basic_002.tcl` source `vmd_view_complexI_9TI4_basic_003.tcl`.
- `vmd_compare_complexI_9TI4_WT_vs_M64V.tcl` sources `vmd_compare_complexI_9TI4_WT_vs_M64V_basic.tcl`.
