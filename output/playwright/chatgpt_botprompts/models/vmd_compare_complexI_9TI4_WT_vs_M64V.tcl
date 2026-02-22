# Wrapper: use the toggle-enabled basic compare view.
#
# Kept for backwards compatibility with older command lines / filenames.
#
# Usage:
#   vmd -e vmd_compare_complexI_9TI4_WT_vs_M64V.tcl -args WT.pdb M64V.pdb --nd ND6 --show nd,lipids

set script_dir [file dirname [info script]]
source [file join $script_dir "vmd_compare_complexI_9TI4_WT_vs_M64V_basic.tcl"]

