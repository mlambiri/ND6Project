# VMD visualization: human Complex I (9TI4) WT view with CLI toggles
#
# Usage (defaults):
#   vmd -e vmd_view_complexI_9TI4_WT.tcl
#
# Usage (explicit file):
#   vmd -e vmd_view_complexI_9TI4_WT.tcl -args /path/to/model.pdb
#
# Usage (options; after -args):
#   --nd ND6 --show context,nd
#   --arm membrane --show context,nd,arm
#   --show all
#
# Interactive toggles (VMD Tk Console):
#   ci_list
#   ci_focus_nd ND6
#   ci_show/ci_hide lipids|cofactors|fes|arm|context|nd

set _ci_script_dir [file dirname [info script]]
source [file join $_ci_script_dir "ci_showhide.tcl"]
unset _ci_script_dir

proc _ci_force_display_prefs {} {
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

proc main {} {
  global argc argv

  set cli_opts [ci_parse_cli $argc $argv]
  if {[dict get $cli_opts help]} {
    puts "Usage: vmd -e vmd_view_complexI_9TI4_WT.tcl -args [pdb] --nd ND6 --show context,nd"
    return
  }

  # Default focus for this script: ND6 unless user specifies --nd.
  if {![dict get $cli_opts nd_explicit]} {
    dict set cli_opts nd_raw ND6
  }
  set vis [ci_resolve_visibility $cli_opts {context nd} {context nd arm lipids cofactors fes}]
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
    set pdb_file [file join $script_dir "complexI_9TI4_WT_heavy_proteinOnly.pdb"]
  }

  set molid [mol new $pdb_file type pdb waitfor all]

  catch {graphics $molid delete all}
  _ci_clear_all_reps $molid
  _ci_force_display_prefs
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

  # Start from a consistent camera and orientation.
  display resetview
  scale by 0.92

  set rot {{0 0 -1 0} {-1 0 0 0} {0 1 0 0} {0 0 0 1}}
  catch {molinfo $molid set rotate_matrix $rot}

  # Context protein.
  ci_add_rep_register "context" $molid {NewCartoon} "protein" {Chain} Opaque

  # ND1..ND6 (tubes). Extra ND6 detail reps are also registered under nd_nd6.
  foreach def [ci_nd_defs] {
    lassign $def name chainid colorid
    set key "nd_[string tolower $name]"
    ci_add_rep_register $key $molid {Tube 0.60 12.0} "protein and chain $chainid" [list ColorID $colorid] Opaque
  }

  # ND6 detailed licorice + residue 64 highlight (WT MET).
  ci_add_rep_register "nd_nd6" $molid {Licorice 0.2 12.0 12.0} "protein and chain m" {ColorID 0} Opaque
  ci_add_rep_register "nd_nd6" $molid {Bonds 0.35 12.0} "protein and chain m and resid 64" {ColorID 4} Opaque

  # Arm representations (membrane-near vs peripheral).
  lassign [ci_make_arm_selections $molid $membrane_dist] sel_arm_mem sel_arm_per
  if {!$arm_include_nd} {
    set sel_arm_mem "$sel_arm_mem and not (chain s i j r l m)"
    set sel_arm_per "$sel_arm_per and not (chain s i j r l m)"
  }
  ci_add_rep_register "arm_membrane" $molid {Tube 0.35 12.0} $sel_arm_mem {ColorID 2} Transparent
  ci_add_rep_register "arm_peripheral" $molid {Tube 0.35 12.0} $sel_arm_per {ColorID 6} Transparent

  # Lipids + cofactors + FeS (only present in full heavy models).
  ci_add_rep_register "lipids" $molid {Bonds 0.20 12.0} "resname CDL PEE PLX DGT" {Resname} Opaque
  ci_add_rep_register "cofactors" $molid {Bonds 0.25 12.0} "resname FMN NDP 8Q1 SF4 FES" {Resname} Opaque
  ci_add_rep_register "fes" $molid {VDW 0.80 12.0} "resname SF4 FES" {Resname} Opaque

  ci_apply_visibility $show_parts $nd_list $arm_mode

  after idle _ci_force_display_prefs
  after 200 _ci_force_display_prefs

  return
}

main
