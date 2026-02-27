# VMD visualization: compare 9TI4 WT vs MT-ND4 R340H (LHON m.11778G>A) (basic style + toggles)
#
# Goal: chain-colored lines with emphasized ND1..ND6 as thick tubes,
# and an ND4 atomistic view + mutation-site highlight for clean WT↔mutant switching.
#
# Usage (defaults):
#   vmd -e vmd_compare_complexI_9TI4_WT_vs_ND4_R340H_basic.tcl
#
# Usage (explicit files):
#   vmd -e vmd_compare_complexI_9TI4_WT_vs_ND4_R340H_basic.tcl -args WT.pdb ND4_R340H.pdb
#
# Usage (options; after -args):
#   --nd ND4 --show nd
#   --show all --arm both
#   --strain wt|mut|both
#
# Interactive toggles (VMD Tk Console):
#   ci_help
#   ci_focus_nd ND4
#   ci_strain wt / ci_strain mut / ci_toggle_strain

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

proc main {} {
  global argc argv

  set cli_opts [ci_parse_cli $argc $argv]
  if {[dict get $cli_opts help]} {
    puts "Usage: vmd -e vmd_compare_complexI_9TI4_WT_vs_ND4_R340H_basic.tcl -args WT.pdb ND4_R340H.pdb --nd ND4 --show nd --strain mut"
    return
  }
  set vis [ci_resolve_visibility $cli_opts]
  set show_parts [dict get $vis show_parts]
  set nd_list [dict get $vis nd_list]
  set arm_mode [dict get $vis arm_mode]
  set arm_include_nd [dict get $vis arm_include_nd]
  set membrane_dist [dict get $vis membrane_dist]
  set strain_mode [dict get $cli_opts strain_mode]
  if {![dict get $cli_opts strain_explicit]} {
    # Default to mutant-only to avoid z-fighting when the two structures overlap.
    set strain_mode mut
  }

  set script_dir [file dirname [info script]]
  set wt_file [dict get $cli_opts wt_file]
  set mut_file [dict get $cli_opts mut_file]
  set positionals [dict get $cli_opts positionals]
  if {$wt_file eq "" && [llength $positionals] >= 1} {
    set wt_file [lindex $positionals 0]
  }
  if {$mut_file eq "" && [llength $positionals] >= 2} {
    set mut_file [lindex $positionals 1]
  }
  if {$wt_file eq "" || $mut_file eq ""} {
    set wt_file  [file join $script_dir "complexI_9TI4_WT_heavy.pdb"]
    set mut_file [file join $script_dir "complexI_9TI4_ND4_R340H_heavy.pdb"]
    if {![file exists $wt_file] || ![file exists $mut_file]} {
      set wt_file  [file join $script_dir "complexI_9TI4_WT_heavy_proteinOnly.pdb"]
      set mut_file [file join $script_dir "complexI_9TI4_ND4_R340H_heavy_proteinOnly.pdb"]
    }
  }

  set wt  [mol new $wt_file  type pdb waitfor all]
  set mut [mol new $mut_file type pdb waitfor all]
  set ::ci_wt_molid $wt
  set ::ci_mut_molid $mut

  catch {graphics $wt delete all}
  catch {graphics $mut delete all}
  _ci_clear_all_reps $wt
  _ci_clear_all_reps $mut

  _ci_force_display_prefs_basic
  ci_reset

  # Align mutant onto WT using CA atoms.
  set sel_wt_fit  [atomselect $wt  "name CA"]
  set sel_mut_fit [atomselect $mut "name CA"]
  if {[$sel_wt_fit num] == [$sel_mut_fit num] && [$sel_wt_fit num] > 0} {
    set M [measure fit $sel_mut_fit $sel_wt_fit]
    set sel_mut_all [atomselect $mut "all"]
    $sel_mut_all move $M
    $sel_mut_all delete
  } else {
    puts "WARNING: Could not align (CA atom counts differ or selection empty)."
  }
  $sel_wt_fit delete
  $sel_mut_fit delete

  # Recenter near origin (based on mutant protein CA).
  set sel_center [atomselect $mut "protein and name CA"]
  if {[$sel_center num] > 0} {
    set pcenter [measure center $sel_center]
    set shift [list \
      [expr {-1.0*[lindex $pcenter 0]}] \
      [expr {-1.0*[lindex $pcenter 1]}] \
      [expr {-1.0*[lindex $pcenter 2]}] \
    ]
    set sel_wt_all  [atomselect $wt  "all"]
    set sel_mut_all [atomselect $mut "all"]
    $sel_wt_all  moveby $shift
    $sel_mut_all moveby $shift
    $sel_wt_all delete
    $sel_mut_all delete
  }
  $sel_center delete

  display resetview
  scale by 0.92

  # Force consistent orientation (same as other scripts).
  set rot {{0 0 -1 0} {-1 0 0 0} {0 1 0 0} {0 0 0 1}}
  catch {molinfo $wt  set rotate_matrix $rot}
  catch {molinfo $mut set rotate_matrix $rot}

  # Context: all protein as chain-colored lines (same style as vmd_view_complexI_9TI4_basic_003.tcl).
  set rep_wt_context [ci_add_rep $wt {Lines} "protein" {Chain} Opaque]
  ci_register_rep "context" $wt $rep_wt_context
  ci_register_rep "wt" $wt $rep_wt_context

  set rep_mut_context [ci_add_rep $mut {Lines} "protein" {Chain} Opaque]
  ci_register_rep "context" $mut $rep_mut_context
  ci_register_rep "mut" $mut $rep_mut_context

  # ND1..ND6 as thicker tubes (per-chain colors like basic_003).
  foreach def [ci_nd_defs] {
    lassign $def name chainid colorid
    set key "nd_[string tolower $name]"

    set rep_wt [ci_add_rep $wt {Tube 0.60 12.0} "protein and chain $chainid" [list ColorID $colorid] Opaque]
    ci_register_rep $key $wt $rep_wt
    ci_register_rep "wt" $wt $rep_wt

    set rep_mut [ci_add_rep $mut {Tube 0.60 12.0} "protein and chain $chainid" [list ColorID $colorid] Opaque]
    ci_register_rep $key $mut $rep_mut
    ci_register_rep "mut" $mut $rep_mut
  }

  # ND4: show *all atoms* (not just the tube) when you do "ci_show ND4".
  # Register the atomistic reps under nd_nd4 so they toggle with ND4.
  set rep_wt_nd4_atoms [ci_add_rep $wt {Licorice 0.20 12.0 12.0} "chain r" {ColorID 7} Opaque]
  ci_register_rep "nd_nd4" $wt $rep_wt_nd4_atoms
  ci_register_rep "wt" $wt $rep_wt_nd4_atoms

  set rep_mut_nd4_atoms [ci_add_rep $mut {Licorice 0.20 12.0 12.0} "chain r" {ColorID 7} Opaque]
  ci_register_rep "nd_nd4" $mut $rep_mut_nd4_atoms
  ci_register_rep "mut" $mut $rep_mut_nd4_atoms

  # Highlight only the *mutation difference* in red on the mutant (new HIS ring atoms).
  set rep_mut_nd4_diff [ci_add_rep $mut {Licorice 0.28 12.0 12.0} "chain r and resid 340 and name ND1 CD2 CE1 NE2" {ColorID 1} Opaque]
  ci_register_rep "nd_nd4" $mut $rep_mut_nd4_diff
  ci_register_rep "mut" $mut $rep_mut_nd4_diff

  # Arm representations (WT + mutant).
  foreach pair [list [list wt $wt] [list mut $mut]] {
    lassign $pair tag molid
    lassign [ci_make_arm_selections $molid $membrane_dist] sel_arm_mem sel_arm_per
    if {!$arm_include_nd} {
      set sel_arm_mem "$sel_arm_mem and not (chain s i j r l m)"
      set sel_arm_per "$sel_arm_per and not (chain s i j r l m)"
    }
    set rep_arm_mem [ci_add_rep $molid {Tube 0.35 12.0} $sel_arm_mem {ColorID 2} Transparent]
    ci_register_rep "arm_membrane" $molid $rep_arm_mem
    ci_register_rep $tag $molid $rep_arm_mem
    set rep_arm_per [ci_add_rep $molid {Tube 0.35 12.0} $sel_arm_per {ColorID 6} Transparent]
    ci_register_rep "arm_peripheral" $molid $rep_arm_per
    ci_register_rep $tag $molid $rep_arm_per
  }

  # Lipids + cofactors/clusters as bonds; Fe-S clusters as tiny spheres (WT + mutant).
  foreach pair [list [list wt $wt] [list mut $mut]] {
    lassign $pair tag molid

    set rep_lipids [ci_add_rep $molid {Bonds 0.20 12.0} "resname CDL PEE PLX DGT" {Resname} Opaque]
    ci_register_rep "lipids" $molid $rep_lipids
    ci_register_rep $tag $molid $rep_lipids

    set rep_cof [ci_add_rep $molid {Bonds 0.25 12.0} "resname FMN NDP 8Q1 SF4 FES" {Resname} Opaque]
    ci_register_rep "cofactors" $molid $rep_cof
    ci_register_rep $tag $molid $rep_cof

    set rep_fes [ci_add_rep $molid {VDW 0.80 12.0} "resname SF4 FES" {Resname} Opaque]
    ci_register_rep "fes" $molid $rep_fes
    ci_register_rep $tag $molid $rep_fes
  }

  ci_apply_visibility $show_parts $nd_list $arm_mode
  ci_strain $strain_mode

  after idle _ci_force_display_prefs_basic
  after 200 _ci_force_display_prefs_basic

  return
}

main
