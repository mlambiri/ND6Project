# Wrapper: use the toggle-enabled basic compare view (WT vs MT-ND1 A52T).
#
# Usage:
#   vmd -e vmd_compare_complexI_9TI4_WT_vs_ND1_A52T.tcl
#   vmd -e vmd_compare_complexI_9TI4_WT_vs_ND1_A52T.tcl -args WT.pdb ND1_A52T.pdb --nd ND1 --show nd --strain mut

set script_dir [file dirname [info script]]
source [file join $script_dir "vmd_compare_complexI_9TI4_WT_vs_ND1_A52T_basic.tcl"]

