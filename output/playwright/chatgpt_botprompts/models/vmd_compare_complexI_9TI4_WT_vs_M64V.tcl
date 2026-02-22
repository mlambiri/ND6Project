# VMD visualization: compare 9TI4 WT vs ND6 M64V (LHON m.14484T>C) models
#
# Usage (defaults):
#   vmd -e vmd_compare_complexI_9TI4_WT_vs_M64V.tcl
#
# Usage (explicit files):
#   vmd -e vmd_compare_complexI_9TI4_WT_vs_M64V.tcl -args /path/WT.pdb /path/M64V.pdb
#
# Visual scheme:
#   - WT protein: gray lines
#   - Mutant protein: chain-colored lines
#   - ND6 chain m: WT blue licorice; mutant red licorice
#   - Residue 64: highlighted (no spheres/VDW) in yellow for both

proc main {} {
  global argc argv

  set script_dir [file dirname [info script]]

  if {$argc >= 2} {
    set wt_file  [lindex $argv 0]
    set mut_file [lindex $argv 1]
  } else {
    # Default to protein-only to avoid noisy PDB plugin warnings about hetero residues.
    set wt_file  [file join $script_dir "complexI_9TI4_WT_heavy_proteinOnly.pdb"]
    set mut_file [file join $script_dir "complexI_9TI4_ND6_M64V_heavy_proteinOnly.pdb"]
    if {![file exists $wt_file] || ![file exists $mut_file]} {
      set wt_file  [file join $script_dir "complexI_9TI4_WT_heavy.pdb"]
      set mut_file [file join $script_dir "complexI_9TI4_ND6_M64V_heavy.pdb"]
    }
  }

  # Load with default bond guessing so licorice/lines work and residues are recognized.
  set wt  [mol new $wt_file  type pdb waitfor all]
  set mut [mol new $mut_file type pdb waitfor all]

  proc _ci_force_display_prefs {} {
    display projection Orthographic
    axes location lowerleft
    display depthcue off
    catch {display backgroundgradient off}
    color Display Background white
    catch {color Display BackgroundTop white}
    catch {color Display BackgroundBottom white}
  }
  _ci_force_display_prefs

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

  # Recenter coordinates near the origin so resetview doesn't throw the complex into a corner.
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

  # Start from a consistent camera and orientation.
  display resetview
  scale by 0.92

  # Force a consistent orientation (axis triad + scene):
  #   - +Z points up on the screen
  #   - +Y points left on the screen
  #   - +X points into the screen (view direction)
  #
  # NOTE: VMD's rotate_matrix convention is the transpose of the intuitive world->view mapping,
  # so we provide the matrix below in VMD's convention.
  set rot {{0 0 -1 0} {-1 0 0 0} {0 1 0 0} {0 0 0 1}}
  catch {molinfo $wt  set rotate_matrix $rot}
  catch {molinfo $mut set rotate_matrix $rot}

  # Clear default reps.
  catch {mol delrep 0 $wt}
  catch {mol delrep 0 $mut}

  # WT: protein (gray ribbons).
  mol representation NewCartoon
  mol selection "protein"
  mol color ColorID 8
  mol material Opaque
  mol addrep $wt

  # WT: ND6 (chain m) highlight (blue).
  mol representation Licorice 0.2 12.0 12.0
  mol selection "chain m"
  mol color ColorID 0
  mol material Opaque
  mol addrep $wt

  # WT: residue 64 highlight (no spheres/VDW).
  mol representation Bonds 0.35 12.0
  mol selection "chain m and resid 64"
  mol color ColorID 4
  mol material Opaque
  mol addrep $wt

  # Mutant: protein (chain-colored ribbons).
  mol representation NewCartoon
  mol selection "protein"
  mol color Chain
  mol material Opaque
  mol addrep $mut

  # Mutant: ND6 (chain m) highlight (red).
  mol representation Licorice 0.2 12.0 12.0
  mol selection "chain m"
  mol color ColorID 1
  mol material Opaque
  mol addrep $mut

  # Mutant: residue 64 highlight (no spheres/VDW).
  mol representation Bonds 0.35 12.0
  mol selection "chain m and resid 64"
  mol color ColorID 4
  mol material Opaque
  mol addrep $mut

  after idle _ci_force_display_prefs
  after 200 _ci_force_display_prefs

  return
}

main
