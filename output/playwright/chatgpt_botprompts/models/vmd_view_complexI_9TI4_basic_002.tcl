# Wrapper: use the toggle-enabled basic view.
#
# Kept for backwards compatibility with older command lines / filenames.
#
# Usage:
#   vmd -e vmd_view_complexI_9TI4_basic_002.tcl -args [pdb] --nd ND6 --arm membrane --show nd,arm,lipids

set script_dir [file dirname [info script]]
source [file join $script_dir "vmd_view_complexI_9TI4_basic_003.tcl"]

