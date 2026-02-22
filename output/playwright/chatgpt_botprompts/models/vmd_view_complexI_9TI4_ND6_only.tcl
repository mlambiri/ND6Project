# VMD visualization: ND6-focused view with CLI toggles (9TI4 exported models)
#
# Usage (defaults to ND6 M64V proteinOnly if present, else WT proteinOnly):
#   vmd -e vmd_view_complexI_9TI4_ND6_only.tcl
#
# Usage (explicit file):
#   vmd -e vmd_view_complexI_9TI4_ND6_only.tcl -args /path/to/model.pdb
#
# Examples:
#   vmd -e vmd_view_complexI_9TI4_ND6_only.tcl -args --show nd --nd ND6
#   vmd -e vmd_view_complexI_9TI4_ND6_only.tcl -args --show nd,lipids --nd ND6
#   vmd -e vmd_view_complexI_9TI4_ND6_only.tcl -args --show all --nd all --arm both
#
# Interactive toggles (VMD Tk Console):
#   ci_list
#   ci_focus_nd ND6
#   ci_show/ci_hide lipids|cofactors|fes|arm|context|nd

set _ci_script_dir [file dirname [info script]]
source [file join $_ci_script_dir "ci_showhide.tcl"]
unset _ci_script_dir

proc _ci_force_display_prefs_nd6 {} {
  display projection Orthographic
  axes location lowerleft
  display depthcue off
  catch {display backgroundgradient off}
  color Display Background white
  catch {color Display BackgroundTop white}
  catch {color Display BackgroundBottom white}
}

proc _ci_clear_all_reps {molid} {
  set n [molinfo $molid get numreps]
  for {set i [expr {$n - 1}]} {$i >= 0} {incr i -1} {
    catch {mol delrep $i $molid}
  }
}

proc _ci_chainids_for_nd_list {nd_list} {
  set chainids {}
  foreach def [ci_nd_defs] {
    lassign $def name chainid colorid
    if {[lsearch -exact $nd_list $name] >= 0} {
      lappend chainids $chainid
    }
  }
  return $chainids
}

proc main {} {
  global argc argv

  set cli_opts [ci_parse_cli $argc $argv]
  if {[dict get $cli_opts help]} {
    puts "Usage: vmd -e vmd_view_complexI_9TI4_ND6_only.tcl -args [pdb] --nd ND6 --show nd,lipids"
    return
  }

  # Default focus for this script: ND6 unless user specifies --nd.
  if {![dict get $cli_opts nd_explicit]} {
    dict set cli_opts nd_raw ND6
  }
  set vis [ci_resolve_visibility $cli_opts {nd} {context nd arm lipids cofactors fes}]
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
    set pdb_file [file join $script_dir "complexI_9TI4_ND6_M64V_heavy_proteinOnly.pdb"]
    if {![file exists $pdb_file]} {
      set pdb_file [file join $script_dir "complexI_9TI4_WT_heavy_proteinOnly.pdb"]
    }
    if {![file exists $pdb_file]} {
      set pdb_file [file join $script_dir "complexI_9TI4_ND6_M64V_heavy.pdb"]
      if {![file exists $pdb_file]} {
        set pdb_file [file join $script_dir "complexI_9TI4_WT_heavy.pdb"]
      }
    }
  }

  set molid [mol new $pdb_file type pdb waitfor all]

  # Start clean if sourced in an existing VMD session.
  catch {graphics $molid delete all}
  _ci_clear_all_reps $molid
  _ci_force_display_prefs_nd6
  ci_reset

  # Center on selected ND subunit(s) to keep the focus in-frame.
  set chainids [_ci_chainids_for_nd_list $nd_list]
  if {[llength $chainids] < 1} {
    set chainids {m}
  }
  set chain_sel [join $chainids " "]
  set sel_center [atomselect $molid "protein and chain $chain_sel and name CA"]
  if {[$sel_center num] < 1} {
    $sel_center delete
    set sel_center [atomselect $molid "protein and chain $chain_sel"]
  }
  if {[$sel_center num] > 0} {
    set c [measure center $sel_center]
    set shift [list \
      [expr {-1.0*[lindex $c 0]}] \
      [expr {-1.0*[lindex $c 1]}] \
      [expr {-1.0*[lindex $c 2]}] \
    ]
    set sel_all [atomselect $molid "all"]
    $sel_all moveby $shift
    $sel_all delete
  } else {
    puts "WARNING: focus selection empty; centering skipped."
  }
  $sel_center delete

  display resetview
  scale by 1.15

  # Match the same orientation convention used in the other scripts:
  #   - +Z up, +Y left, +X into screen.
  set rot {{0 0 -1 0} {-1 0 0 0} {0 1 0 0} {0 0 0 1}}
  catch {molinfo $molid set rotate_matrix $rot}

  # Context: all protein as chain-colored lines (off by default).
  ci_add_rep_register "context" $molid {Lines} "protein" {Chain} Opaque

  # ND1..ND6 tubes (so you can switch focus with --nd / ci_focus_nd).
  foreach def [ci_nd_defs] {
    lassign $def name chainid colorid
    set key "nd_[string tolower $name]"
    ci_add_rep_register $key $molid {Tube 0.60 12.0} "protein and chain $chainid" [list ColorID $colorid] Opaque
  }

  # ND6 extra detail reps (thick tube + licorice + resid 64 highlight).
  ci_add_rep_register "nd_nd6" $molid {Tube 0.90 12.0} "protein and chain m" {Chain} Opaque
  ci_add_rep_register "nd_nd6" $molid {Licorice 0.22 12.0 12.0} "protein and chain m" {Chain} Opaque
  ci_add_rep_register "nd_nd6" $molid {Bonds 0.45 12.0} "protein and chain m and resid 64" {ColorID 4} Opaque

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

  after idle _ci_force_display_prefs_nd6
  after 200 _ci_force_display_prefs_nd6

  return
}

main
