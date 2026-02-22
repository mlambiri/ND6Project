# Wrapper: use the toggle-enabled compare-features view.
#
# Kept for backwards compatibility with older command lines / filenames.
#
# Usage:
#   vmd -e vmd_pro.tcl -args WT.pdb M64V.pdb --show all --arm both

set script_dir [file dirname [info script]]
source [file join $script_dir "vmd_compare_complexI_9TI4_WT_vs_M64V_features.tcl"]

