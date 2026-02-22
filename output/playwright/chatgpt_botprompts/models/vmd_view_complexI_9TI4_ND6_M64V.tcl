# VMD visualization: human Complex I (9TI4) ND6 M64V heavy-atoms PDB
#
# Usage:
#   vmd -e vmd_view_complexI_9TI4_ND6_M64V.tcl
#
# Notes:
#   - Expects to be run from this folder (or it resolves paths relative to the script).
#   - Highlights MT-ND6 as chain `m` and residue 64 (mutant VAL) in yellow.

proc main {} {
  global argc argv

  set script_dir [file dirname [info script]]
  if {$argc >= 1} {
    set pdb_file [lindex $argv 0]
  } else {
    set pdb_file [file join $script_dir "complexI_9TI4_ND6_M64V_heavy_proteinOnly.pdb"]
  }

  mol new $pdb_file type pdb waitfor all
  set molid [molinfo top]

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

  # Force a consistent orientation (axis triad + scene):
  #   - +Z points up on the screen
  #   - +Y points left on the screen
  #   - +X points into the screen (view direction)
  #
  # NOTE: VMD's rotate_matrix convention is the transpose of the intuitive world->view mapping,
  # so we provide the matrix below in VMD's convention.
  set rot {{0 0 -1 0} {-1 0 0 0} {0 1 0 0} {0 0 0 1}}
  catch {molinfo $molid set rotate_matrix $rot}

  # Remove default representation.
  catch {mol delrep 0 $molid}

  # Complex I protein (ribbons).
  mol representation NewCartoon
  mol selection "protein"
  mol color Chain
  mol material Opaque
  mol addrep $molid

  # ND6 (chain m).
  mol representation Licorice 0.2 12.0 12.0
  mol selection "chain m"
  mol color ColorID 1
  mol material Opaque
  mol addrep $molid

  # ND6 residue 64 (mutant VAL).
  mol representation Bonds 0.35 12.0
  mol selection "chain m and resid 64"
  mol color ColorID 4
  mol material Opaque
  mol addrep $molid

  # Optional label at CA.
  set sel_ca [atomselect $molid "chain m and resid 64 and name CA"]
  if {[$sel_ca num] == 1} {
    label add Atoms $molid [lindex [$sel_ca get index] 0]
  }
  $sel_ca delete

  after idle _ci_force_display_prefs
  after 200 _ci_force_display_prefs

  return
}

main
