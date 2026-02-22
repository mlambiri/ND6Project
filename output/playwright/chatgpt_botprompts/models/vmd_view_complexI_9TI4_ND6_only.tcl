# VMD visualization: show ONLY MT-ND6 (9TI4 exported models: chain m)
#
# Usage (defaults to ND6 M64V proteinOnly if present, else WT proteinOnly):
#   vmd -e vmd_view_complexI_9TI4_ND6_only.tcl
#
# Usage (explicit file):
#   vmd -e vmd_view_complexI_9TI4_ND6_only.tcl -args /path/to/model.pdb
#
# Notes:
# - ND6 is chain `m` in the PDBs generated in this folder.
# - Shows ND6 as a thick tube + licorice, and highlights resid 64.

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

proc main {} {
  global argc argv

  set script_dir [file dirname [info script]]
  if {$argc >= 1} {
    set pdb_file [lindex $argv 0]
  } else {
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

  # Center on ND6 to keep it in-frame.
  set sel_center [atomselect $molid "protein and chain m and name CA"]
  if {[$sel_center num] < 1} {
    $sel_center delete
    set sel_center [atomselect $molid "protein and chain m"]
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
    puts "WARNING: ND6 selection empty; centering skipped. (Is ND6 still chain m in this PDB?)"
  }
  $sel_center delete

  display resetview
  scale by 1.15

  # Match the same orientation convention used in the other scripts:
  #   - +Z up, +Y left, +X into screen.
  set rot {{0 0 -1 0} {-1 0 0 0} {0 1 0 0} {0 0 0 1}}
  catch {molinfo $molid set rotate_matrix $rot}

  # ND6 backbone as a thick tube.
  mol representation Tube 0.90 12.0
  mol selection "protein and chain m"
  mol color Chain
  mol material Opaque
  mol addrep $molid

  # ND6 atoms (for context).
  mol representation Licorice 0.22 12.0 12.0
  mol selection "protein and chain m"
  mol color Chain
  mol material Opaque
  mol addrep $molid

  # Highlight residue 64 (WT MET / mutant VAL) with bonds.
  mol representation Bonds 0.45 12.0
  mol selection "protein and chain m and resid 64"
  mol color ColorID 4
  mol material Opaque
  mol addrep $molid

  after idle _ci_force_display_prefs_nd6
  after 200 _ci_force_display_prefs_nd6

  return
}

main

