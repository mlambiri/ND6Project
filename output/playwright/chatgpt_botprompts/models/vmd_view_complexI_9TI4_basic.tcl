# VMD visualization: human Complex I (9TI4) basic view (no labels/callouts)
#
# Goal: resemble the "chain-colored lines + highlighted membrane arm" style,
# but WITHOUT any big spheres and WITHOUT any text.
#
# Usage (defaults to ND6 M64V heavy model if present, else WT):
#   vmd -e vmd_view_complexI_9TI4_basic.tcl
#
# Usage (explicit file):
#   vmd -e vmd_view_complexI_9TI4_basic.tcl -args /path/to/model.pdb
#
# Notes:
# - This script does not draw any labels/callouts.
# - It avoids VDW/CPK representations (no large spheres).

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
  if {$argc >= 1} {
    set pdb_file [lindex $argv 0]
  } else {
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

  set molid [mol new $pdb_file type pdb waitfor all]

  # Ensure we start clean if sourced in an existing session.
  catch {graphics $molid delete all}
  _ci_clear_all_reps $molid

  _ci_force_display_prefs_basic

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

  # Base: all protein as chain-colored lines.
  mol representation Lines
  mol selection "protein"
  mol color Chain
  mol material Opaque
  mol addrep $molid

  # Emphasize mitochondrial-encoded membrane subunits (9TI4 chain IDs in our exported PDBs):
  #   ND1=s ND2=i ND3=j ND4=r ND5=l ND6=m
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
    mol addrep $molid
  }

  # Show native lipids + cofactors/clusters as bonds/lines (small, no spheres).
  mol representation Bonds 0.20 12.0
  mol selection "resname CDL PEE PLX"
  mol color Resname
  mol material Opaque
  mol addrep $molid

  mol representation Bonds 0.25 12.0
  mol selection "resname FMN NDP 8Q1 SF4 FES"
  mol color Resname
  mol material Opaque
  mol addrep $molid

  # Fe-S clusters as spheres (small selection; avoids the "big spheres" issue).
  mol representation VDW 0.80 12.0
  mol selection "resname SF4 FES"
  mol color Resname
  mol material Opaque
  mol addrep $molid

  # Mark ND6 location with a leader line + text (chain m).
  set sel_nd6 [atomselect $molid "protein and chain m and name CA"]
  if {[$sel_nd6 num] > 0} {
    set c [measure center $sel_nd6]
    $sel_nd6 delete

    set sel_prot [atomselect $molid "protein and name CA"]
    set mm [measure minmax $sel_prot]
    $sel_prot delete
    set pmin [lindex $mm 0]
    set pmax [lindex $mm 1]
    set zmin [lindex $pmin 2]
    set zmax [lindex $pmax 2]
    set zspan [expr {$zmax - $zmin}]
    set zmargin [expr {0.06*$zspan + 12.0}]

    # Place label under the complex, vertically aligned with ND6 (same Y).
    set lpos [list [lindex $c 0] [lindex $c 1] [expr {$zmin - $zmargin}]]
    graphics $molid color black
    graphics $molid line $c $lpos width 2
    _ci_graphics_text $molid $lpos "ND6" 1.6
  } else {
    $sel_nd6 delete
    puts "NOTE: ND6 label skipped (chain m not found)."
  }

  # Re-apply prefs after startup scripts/themes.
  after idle _ci_force_display_prefs_basic
  after 200 _ci_force_display_prefs_basic

  return
}

main
