# Complex I (human) + ND6 M64V starter models

This folder contains "clean" (heavy-atoms-only) coordinate files exported from RCSB **PDB 9TI4** (human mitochondrial Complex I; cryo-EM) and a simple ND6 mutation helper.

## Files

- `complexI_9TI4_WT_heavy.pdb`: full Complex I, heavy atoms only (HOH removed)
- `complexI_9TI4_ND6_M64V_heavy.pdb`: same, but **MT-ND6 MET64→VAL** applied
- `complexI_9TI4_WT_heavy_proteinOnly.pdb`: protein-only (ATOM records), heavy atoms only
- `complexI_9TI4_ND6_M64V_heavy_proteinOnly.pdb`: protein-only (ATOM records), heavy atoms only, mutated
- `nd6_chain_m_WT_heavy.pdb`: ND6 only (chain `m`)
- `nd6_chain_m_M64V_heavy.pdb`: ND6 only (chain `m`), mutated
- `complexI_9TI4_chain_map.json`: author chain IDs → PDB chain IDs (identity for 9TI4)

## Variant mapping

LHON mtDNA variant **m.14484T>C** in **MT-ND6** maps to protein change **p.Met64Val (M64V)**.

In **9TI4**, ND6 corresponds to **chain `m`**, and the mutated residue is **resid 64**.

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

## VMD visualization

From this folder:

```bash
vmd -e vmd_view_complexI_9TI4_WT.tcl
vmd -e vmd_view_complexI_9TI4_ND6_M64V.tcl
vmd -e vmd_compare_complexI_9TI4_WT_vs_M64V.tcl
vmd -e vmd_view_complexI_9TI4_features.tcl
vmd -e vmd_compare_complexI_9TI4_WT_vs_M64V_features.tcl
```

By default, the basic scripts (`vmd_view_*` / `vmd_compare_*`) load the `*_proteinOnly.pdb` files (less noisy in VMD). The `*_features.tcl` scripts default to the full heavy-atom models so lipids/cofactors are visible. Pass explicit PDB paths via `-args` as needed.
