# VMD visualization: compare 9TI4 WT vs ND6 M64V (LHON m.14484T>C) (basic style)
#
# Goal: chain-colored lines with emphasized ND1..ND6 as thick tubes,
# WITHOUT any labels/callouts and WITHOUT any big spheres.
#
# Usage (defaults):
#   vmd -e vmd_compare_complexI_9TI4_WT_vs_M64V_basic.tcl
#
# Usage (explicit files):
#   vmd -e vmd_compare_complexI_9TI4_WT_vs_M64V_basic.tcl -args WT.pdb M64V.pdb

proc _ci_force_display_prefs_basic {} {
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

proc _ci_graphics_text {molid xyz text {size 1.4}} {
  # Black text with a subtle white halo for readability on busy geometry.
  graphics $molid color white
  graphics $molid text $xyz $text size [expr {$size + 0.25}]
  graphics $molid color black
  graphics $molid text $xyz $text size $size
}

proc main {} {
  global argc argv

  set script_dir [file dirname [info script]]
  if {$argc >= 2} {
    set wt_file  [lindex $argv 0]
    set mut_file [lindex $argv 1]
  } else {
    set wt_file  [file join $script_dir "complexI_9TI4_WT_heavy.pdb"]
    set mut_file [file join $script_dir "complexI_9TI4_ND6_M64V_heavy.pdb"]
    if {![file exists $wt_file] || ![file exists $mut_file]} {
      set wt_file  [file join $script_dir "complexI_9TI4_WT_heavy_proteinOnly.pdb"]
      set mut_file [file join $script_dir "complexI_9TI4_ND6_M64V_heavy_proteinOnly.pdb"]
    }
  }

  set wt  [mol new $wt_file  type pdb waitfor all]
  set mut [mol new $mut_file type pdb waitfor all]

  catch {graphics $wt delete all}
  catch {graphics $mut delete all}
  _ci_clear_all_reps $wt
  _ci_clear_all_reps $mut

  _ci_force_display_prefs_basic

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

  # WT: protein as gray lines.
  mol representation Lines
  mol selection "protein"
  mol color ColorID 8
  mol material Opaque
  mol addrep $wt

  # Mutant: protein as chain-colored lines.
  mol representation Lines
  mol selection "protein"
  mol color Chain
  mol material Opaque
  mol addrep $mut

  # Emphasize ND1..ND6 on mutant as thick tubes (no labels).
  set nd_defs {
    {ND1 s 1}
    {ND2 i 3}
    {ND3 j 4}
    {ND4 r 7}
    {ND5 l 0}
    {ND6 m 9}
  }
  foreach def $nd_defs {
    lassign $def name chainid colorid
    mol representation Tube 0.60 12.0
    mol selection "protein and chain $chainid"
    mol color ColorID $colorid
    mol material Opaque
    mol addrep $mut
  }

  # Mutant: lipids + cofactors/clusters as bonds (small, no spheres).
  mol representation Bonds 0.20 12.0
  mol selection "resname CDL PEE PLX"
  mol color Resname
  mol material Opaque
  mol addrep $mut

  mol representation Bonds 0.25 12.0
  mol selection "resname FMN NDP 8Q1 SF4 FES"
  mol color Resname
  mol material Opaque
  mol addrep $mut

  # Fe-S clusters as spheres (small selection; avoids the "big spheres" issue).
  mol representation VDW 0.80 12.0
  mol selection "resname SF4 FES"
  mol color Resname
  mol material Opaque
  mol addrep $mut

  # Mark ND6 location on mutant with a leader line + text (chain m).
  set sel_nd6 [atomselect $mut "protein and chain m and name CA"]
  if {[$sel_nd6 num] > 0} {
    set c [measure center $sel_nd6]
    $sel_nd6 delete

    set sel_prot [atomselect $mut "protein and name CA"]
    set mm [measure minmax $sel_prot]
    $sel_prot delete
    set pmin [lindex $mm 0]
    set pmax [lindex $mm 1]
    set zmin [lindex $pmin 2]
    set zmax [lindex $pmax 2]
    set zspan [expr {$zmax - $zmin}]
    set zmargin [expr {0.06*$zspan + 12.0}]

    set lpos [list [lindex $c 0] [lindex $c 1] [expr {$zmin - $zmargin}]]
    graphics $mut color black
    graphics $mut line $c $lpos width 2
    _ci_graphics_text $mut $lpos "ND6" 1.6
  } else {
    $sel_nd6 delete
    puts "NOTE: ND6 label skipped (chain m not found)."
  }

  after idle _ci_force_display_prefs_basic
  after 200 _ci_force_display_prefs_basic

  return
}

main
