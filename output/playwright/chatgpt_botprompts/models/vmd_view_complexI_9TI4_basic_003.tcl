# VMD visualization: human Complex I (9TI4) basic view with CLI toggles (no labels/callouts)
#
# Goal: resemble the "chain-colored lines + highlighted membrane arm" style,
# but WITHOUT any big spheres and WITHOUT any text.
#
# Usage (defaults to ND6 M64V heavy model if present, else WT):
#   vmd -e vmd_view_complexI_9TI4_basic_003.tcl
#
# Usage (explicit file):
#   vmd -e vmd_view_complexI_9TI4_basic_003.tcl -args /path/to/model.pdb
#
# Usage (options):
#   vmd -e vmd_view_complexI_9TI4_basic_003.tcl -args [pdb] --nd ND6 --arm membrane --show context,nd,arm,lipids
#   vmd -e vmd_view_complexI_9TI4_basic_003.tcl -args --show nd,arm,lipids --nd ND3 --hide cofactors,fes
#
# Options (after -args):
#   [pdb_file]                 Optional first positional arg.
#   --pdb <file.pdb>           Explicit PDB path.
#   --nd <ND1,ND2,...|all|none> ND subunits to highlight (default: all).
#   --show <default|all|csv>   Parts to show: context,nd,arm,lipids,cofactors,fes
#                              (default: context,nd,lipids,cofactors,fes).
#   --hide <csv>               Parts to hide (removes from --show).
#   --arm <membrane|peripheral|both|off>
#                              If set (and --show not provided), also enables 'arm'.
#   --arm-include-nd <0|1>     Include ND chains in the 'arm' representation (default: 0).
#   --membrane-dist <A>        Distance cutoff for membrane arm selection (default: 5.0).
#   --help                     Print help and return.
#
# Interactive toggles (after load, in VMD Tk Console):
#   ci_list
#   ci_show lipids
#   ci_hide lipids
#   ci_focus_nd ND6
#
# Notes:
# - Lipids/cofactors require the full heavy-atom PDB; proteinOnly models won't show them.
# - No labels/callouts are drawn.

set _ci_script_dir [file dirname [info script]]
source [file join $_ci_script_dir "ci_showhide.tcl"]
unset _ci_script_dir

proc _ci_force_display_prefs_basic {} {
  display projection Orthographic
  axes location lowerleft
  display depthcue off
  catch {display backgroundgradient off}
  color Display Background black
  catch {color Display BackgroundTop black}
  catch {color Display BackgroundBottom black}
}

proc _ci_clear_all_reps {molid} {
  set n [molinfo $molid get numreps]
  for {set i [expr {$n - 1}]} {$i >= 0} {incr i -1} {
    catch {mol delrep $i $molid}
  }
}

proc _ci_usage {} {
  puts ""
  puts "vmd_view_complexI_9TI4_basic_003.tcl"
  puts "  vmd -e vmd_view_complexI_9TI4_basic_003.tcl"
  puts "  vmd -e vmd_view_complexI_9TI4_basic_003.tcl -args /path/to/model.pdb"
  puts "  vmd -e vmd_view_complexI_9TI4_basic_003.tcl -args --nd ND6 --arm membrane --show context,nd,arm,lipids"
  puts ""
  puts "Options:"
  puts "  --pdb <file.pdb>"
  puts "  --nd <ND1,ND2,...|all|none>"
  puts "  --show <default|all|csv parts>"
  puts "  --hide <csv parts>"
  puts "  --arm <membrane|peripheral|both|off>"
  puts "  --arm-include-nd <0|1>"
  puts "  --membrane-dist <A>"
  puts "  --help"
  puts ""
  puts "Parts: context, nd, arm, lipids, cofactors, fes"
  puts ""
}

proc main {} {
  global argc argv

  set cli_opts [ci_parse_cli $argc $argv]
  if {[dict get $cli_opts help]} {
    _ci_usage
    return
  }
  set vis [ci_resolve_visibility $cli_opts]
  set show_parts [dict get $vis show_parts]
  set nd_list [dict get $vis nd_list]
  set arm_mode [dict get $vis arm_mode]
  set arm_include_nd [dict get $vis arm_include_nd]
  set membrane_dist [dict get $vis membrane_dist]

  set script_dir [file dirname [info script]]
  set pdb_file [dict get $cli_opts pdb_file]
  set positionals [dict get $cli_opts positionals]
  if {$pdb_file eq "" && [llength $positionals] >= 1} {
    set pdb_file [lindex $positionals 0]
  }
  if {$pdb_file eq ""} {
    # Prefer mutant model if it exists.
    set pdb_file [file join $script_dir "complexI_9TI4_ND6_M64V_heavy.pdb"]
    if {![file exists $pdb_file]} {
      set pdb_file [file join $script_dir "complexI_9TI4_WT_heavy.pdb"]
    }
    if {![file exists $pdb_file]} {
      set pdb_file [file join $script_dir "complexI_9TI4_ND6_M64V_heavy_proteinOnly.pdb"]
      if {![file exists $pdb_file]} {
        set pdb_file [file join $script_dir "complexI_9TI4_WT_heavy_proteinOnly.pdb"]
      }
    }
  }

  puts "Complex I view (basic_003):"
  puts "  PDB: $pdb_file"
  puts "  show parts: [join $show_parts {, }]"
  if {[llength $nd_list] > 0} {
    puts "  ND highlight: [join $nd_list {, }]"
  } else {
    puts "  ND highlight: (none)"
  }
  if {[lsearch -exact $show_parts "arm"] >= 0} {
    puts "  arm: $arm_mode (membrane_dist=$membrane_dist A, include_nd=$arm_include_nd)"
  }

  set molid [mol new $pdb_file type pdb waitfor all]

  # Ensure we start clean if sourced in an existing session.
  catch {graphics $molid delete all}
  _ci_clear_all_reps $molid

  _ci_force_display_prefs_basic

  # Reset shared show/hide registry for this scene.
  ci_reset

  # Recenter coordinates near the origin so resetview doesn't throw the complex into a corner.
  set sel_center [atomselect $molid "protein and name CA"]
  if {[$sel_center num] > 0} {
    set pcenter [measure center $sel_center]
    set shift [list \
      [expr {-1.0*[lindex $pcenter 0]}] \
      [expr {-1.0*[lindex $pcenter 1]}] \
      [expr {-1.0*[lindex $pcenter 2]}] \
    ]
    set sel_all [atomselect $molid "all"]
    $sel_all moveby $shift
    $sel_all delete
  }
  $sel_center delete

  display resetview
  scale by 0.92

  # Force a consistent orientation:
  #   - +Z points up on the screen
  #   - +Y points left on the screen
  #   - +X points into the screen (view direction)
  set rot {{0 0 -1 0} {-1 0 0 0} {0 1 0 0} {0 0 0 1}}
  catch {molinfo $molid set rotate_matrix $rot}

  # --- Create representations (register all; visibility is set afterwards) ---

  # Context: all protein as chain-colored lines.
  ci_add_rep_register "context" $molid {Lines} "protein" {Chain} Opaque

  # ND subunits as thicker tubes (per-chain so colors stay distinct).
  foreach def [ci_nd_defs] {
    lassign $def name chainid colorid
    set key "nd_[string tolower $name]"
    ci_add_rep_register $key $molid {Tube 0.60 12.0} "protein and chain $chainid" [list ColorID $colorid] Opaque
  }

  # Arm representations: membrane-near vs peripheral, using proximity to lipids (or ND neighborhood fallback).
  lassign [ci_make_arm_selections $molid $membrane_dist] sel_arm_mem sel_arm_per
  if {!$arm_include_nd} {
    set sel_arm_mem "$sel_arm_mem and not (chain s i j r l m)"
    set sel_arm_per "$sel_arm_per and not (chain s i j r l m)"
  }
  ci_add_rep_register "arm_membrane" $molid {Tube 0.35 12.0} $sel_arm_mem {ColorID 2} Transparent
  ci_add_rep_register "arm_peripheral" $molid {Tube 0.35 12.0} $sel_arm_per {ColorID 6} Transparent

  # Lipids.
  ci_add_rep_register "lipids" $molid {Bonds 0.20 12.0} "resname CDL PEE PLX DGT" {Resname} Opaque

  # Cofactors/clusters as bonds (small, no spheres).
  ci_add_rep_register "cofactors" $molid {Bonds 0.25 12.0} "resname FMN NDP 8Q1 SF4 FES" {Resname} Opaque

  # Fe-S clusters as spheres (small selection; avoids the "big spheres" issue).
  ci_add_rep_register "fes" $molid {VDW 0.80 12.0} "resname SF4 FES" {Resname} Opaque

  # Apply initial visibility based on CLI toggles.
  ci_apply_visibility $show_parts $nd_list $arm_mode

  # Re-apply prefs after startup scripts/themes.
  after idle _ci_force_display_prefs_basic
  after 200 _ci_force_display_prefs_basic

  return
}

main
