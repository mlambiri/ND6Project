# Wrapper: use the toggle-enabled basic compare view (WT vs MT-ND4 R340H).
#
# Usage:
#   vmd -e vmd_compare_complexI_9TI4_WT_vs_ND4_R340H.tcl
#   vmd -e vmd_compare_complexI_9TI4_WT_vs_ND4_R340H.tcl -args WT.pdb ND4_R340H.pdb --nd ND4 --show nd --strain mut

set script_dir [file dirname [info script]]
source [file join $script_dir "vmd_compare_complexI_9TI4_WT_vs_ND4_R340H_basic.tcl"]

