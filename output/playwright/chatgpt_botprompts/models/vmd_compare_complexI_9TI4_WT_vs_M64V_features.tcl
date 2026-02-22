# VMD visualization: compare 9TI4 WT vs ND6 M64V (LHON m.14484T>C) with key features marked
#
# Usage (defaults):
#   vmd -e vmd_compare_complexI_9TI4_WT_vs_M64V_features.tcl
#
# Usage (explicit files):
#   vmd -e vmd_compare_complexI_9TI4_WT_vs_M64V_features.tcl -args WT.pdb M64V.pdb
#
# Defaults load the *full* heavy-atom models (includes lipids/cofactors):
#   - complexI_9TI4_WT_heavy.pdb
#   - complexI_9TI4_ND6_M64V_heavy.pdb
#
# Feature mapping in 9TI4 (auth chain IDs in our exported PDBs):
#   ND1  = chain s
#   ND2  = chain i
#   ND3  = chain j
#   ND4  = chain r
#   ND5  = chain l
#   ND6  = chain m   (mutation site resid 64)
#
# Fe-S clusters / cofactors (by residue identity present in 9TI4):
#   FMN: resname FMN (chain C resid 502)
#   NADPH: resname NDP (chain J resid 401)
#   Fe-S: resname FES (2Fe-2S), SF4 (4Fe-4S)
#   Quinone-like ligand: resname 8Q1 (present on chains F and X; Q-site likely chain F)
#   Lipids: resname CDL, PEE, PLX
#
# Note: labels for N1a..N6b are assigned based on subunit association and residue IDs in 9TI4.

proc _center_of_selection {molid seltext} {
  set sel [atomselect $molid $seltext]
  if {[$sel num] < 1} {
    $sel delete
    return ""
  }
  set c [measure center $sel]
  $sel delete
  return $c
}

proc _vadd {a b} {
  return [list \
    [expr {[lindex $a 0] + [lindex $b 0]}] \
    [expr {[lindex $a 1] + [lindex $b 1]}] \
    [expr {[lindex $a 2] + [lindex $b 2]}] \
  ]
}

proc _vsub {a b} {
  return [list \
    [expr {[lindex $a 0] - [lindex $b 0]}] \
    [expr {[lindex $a 1] - [lindex $b 1]}] \
    [expr {[lindex $a 2] - [lindex $b 2]}] \
  ]
}

proc _vdot {a b} {
  return [expr {([lindex $a 0]*[lindex $b 0]) + ([lindex $a 1]*[lindex $b 1]) + ([lindex $a 2]*[lindex $b 2])}]
}

proc _mat3_from_rotate {rot} {
  # Accept either {{..4..} x4} or flat 16-element lists; return 3x3 rows.
  if {[llength $rot] == 4 && [llength [lindex $rot 0]] == 4} {
    set r0 [lrange [lindex $rot 0] 0 2]
    set r1 [lrange [lindex $rot 1] 0 2]
    set r2 [lrange [lindex $rot 2] 0 2]
    return [list $r0 $r1 $r2]
  }
  if {[llength $rot] == 16} {
    set r0 [list [lindex $rot 0] [lindex $rot 1] [lindex $rot 2]]
    set r1 [list [lindex $rot 4] [lindex $rot 5] [lindex $rot 6]]
    set r2 [list [lindex $rot 8] [lindex $rot 9] [lindex $rot 10]]
    return [list $r0 $r1 $r2]
  }
  # Fallback: assume 3 rows of 3.
  return [list [lindex $rot 0] [lindex $rot 1] [lindex $rot 2]]
}

proc _mol_to_world {rows v} {
  set r0 [lindex $rows 0]
  set r1 [lindex $rows 1]
  set r2 [lindex $rows 2]
  return [list [_vdot $r0 $v] [_vdot $r1 $v] [_vdot $r2 $v]]
}

proc _world_to_mol {rows w} {
  # For a pure rotation R, inverse is transpose: mol = R^T * world.
  lassign $w x y z
  set r0 [lindex $rows 0]
  set r1 [lindex $rows 1]
  set r2 [lindex $rows 2]
  set mx [expr {$x*[lindex $r0 0] + $y*[lindex $r1 0] + $z*[lindex $r2 0]}]
  set my [expr {$x*[lindex $r0 1] + $y*[lindex $r1 1] + $z*[lindex $r2 1]}]
  set mz [expr {$x*[lindex $r0 2] + $y*[lindex $r1 2] + $z*[lindex $r2 2]}]
  return [list $mx $my $mz]
}

proc _world_extents_from_minmax {rows center minv maxv} {
  # Transform the 8 corners of an axis-aligned box to world coords; return min/max per axis.
  set wxmin 1e99
  set wymin 1e99
  set wzmin 1e99
  set wxmax -1e99
  set wymax -1e99
  set wzmax -1e99
  foreach x [list [lindex $minv 0] [lindex $maxv 0]] {
    foreach y [list [lindex $minv 1] [lindex $maxv 1]] {
      foreach z [list [lindex $minv 2] [lindex $maxv 2]] {
        set p [list $x $y $z]
        set w [_mol_to_world $rows [_vsub $p $center]]
        set wx [lindex $w 0]
        set wy [lindex $w 1]
        set wz [lindex $w 2]
        if {$wx < $wxmin} { set wxmin $wx }
        if {$wy < $wymin} { set wymin $wy }
        if {$wz < $wzmin} { set wzmin $wz }
        if {$wx > $wxmax} { set wxmax $wx }
        if {$wy > $wymax} { set wymax $wy }
        if {$wz > $wzmax} { set wzmax $wz }
      }
    }
  }
  return [list $wxmin $wxmax $wymin $wymax $wzmin $wzmax]
}

proc _graphics_text {molid xyz text colorname {size 2.5}} {
  # Readable label with a simple white outline (masks busy geometry behind text).
  if {$xyz eq ""} { return }
  graphics $molid color white
  graphics $molid text $xyz $text size [expr {$size + 0.20}]
  graphics $molid color $colorname
  graphics $molid text $xyz $text size $size
}

proc _label_at_atom {molid seltext label colorname} {
  set sel [atomselect $molid $seltext]
  if {[$sel num] < 1} {
    puts "NOTE: missing label target: $seltext"
    $sel delete
    return
  }
  set xyz [lindex [$sel get {x y z}] 0]
  $sel delete
  _graphics_text $molid $xyz $label $colorname 0.9
}

proc _set_axis_component {xyz axis value} {
  set out $xyz
  lset out $axis $value
  return $out
}

proc _pack_positions {pairs minsep} {
  # Given a list of {key desiredPos} pairs, return {key packedPos} with minimum spacing.
  set sorted [lsort -real -index 1 $pairs]
  set packed {}
  set prev -1e99
  foreach p $sorted {
    lassign $p key pos
    if {$pos < ($prev + $minsep)} {
      set pos [expr {$prev + $minsep}]
    }
    lappend packed [list $key $pos]
    set prev $pos
  }
  return $packed
}

proc _callout {molid from_xyz to_xyz text colorname size} {
  if {$from_xyz eq "" || $to_xyz eq ""} { return }
  # Callout with the *label-end* segment aligned to the Cartesian Y axis
  # (in this view: Z up, Y left/right, X into the screen).
  set elbow [list [lindex $to_xyz 0] [lindex $from_xyz 1] [lindex $to_xyz 2]]
  graphics $molid color black
  graphics $molid line $from_xyz $elbow width 2
  graphics $molid line $elbow $to_xyz width 2
  _graphics_text $molid $to_xyz $text $colorname $size
}

proc _ci_force_display_prefs {} {
  # Some user VMD startup scripts can override display prefs; enforce ours again.
  display projection Orthographic
  axes location lowerleft
  display depthcue off
  catch {display backgroundgradient off}
  color Display Background white
  catch {color Display BackgroundTop white}
  catch {color Display BackgroundBottom white}
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

  # Clear any lingering OpenGL drawing objects (helps if you source the script in an existing VMD session).
  catch {graphics $wt delete all}
  catch {graphics $mut delete all}

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
  set pcenter [_center_of_selection $mut "protein and name CA"]
  if {$pcenter ne ""} {
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

  # Start from a consistent camera.
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
  if {[catch {molinfo $wt  set rotate_matrix $rot} err]} {
    puts "WARNING: could not set rotate_matrix on WT: $err"
  }
  if {[catch {molinfo $mut set rotate_matrix $rot} err]} {
    puts "WARNING: could not set rotate_matrix on MUT: $err"
  }

  # Clear default reps.
  catch {mol delrep 0 $wt}
  catch {mol delrep 0 $mut}

  # --- Base representations ---
  # WT: protein in gray (reference).
  mol representation NewCartoon
  mol selection "protein"
  mol color ColorID 8
  mol material Opaque
  mol addrep $wt

  # MUT: protein chain-colored.
  mol representation NewCartoon
  mol selection "protein"
  mol color Chain
  mol material Opaque
  mol addrep $mut

  # --- ND1..ND6 marked on mutant ---
  # ND labels: below the complex (low Z)
  # Feature labels: above the complex (high Z)
  # All label anchors packed along Cartesian Y (horizontal in this view).
  set xaxis 0
  set yaxis 1
  set zaxis 2

  # Use the full protein box for outside placement (below/above).
  set sel_prot [atomselect $mut "protein and name CA"]
  set mm_prot [measure minmax $sel_prot]
  $sel_prot delete
  set pmin [lindex $mm_prot 0]
  set pmax [lindex $mm_prot 1]
  set xmin [lindex $pmin $xaxis]
  set xmax [lindex $pmax $xaxis]
  set zmin [lindex $pmin $zaxis]
  set zmax [lindex $pmax $zaxis]
  set xspan [expr {$xmax - $xmin}]
  set zspan [expr {$zmax - $zmin}]

  # X is "into screen" in this view; put labels slightly toward the viewer (more negative X).
  set xmargin [expr {0.08*$xspan + 6.0}]
  # Keep Z margins modest so labels stay in-frame without extreme zooming.
  set zmargin [expr {0.03*$zspan + 6.0}]

  set label_x [expr {$xmin - $xmargin}]
  set nd_z    [expr {$zmin - $zmargin}]
  set feat_z  [expr {$zmax + $zmargin}]

  # name chain colorid
  set nd_defs {
    {ND1 s 1}
    {ND2 i 3}
    {ND3 j 4}
    {ND4 r 7}
    {ND5 l 0}
    {ND6 m 9}
  }

  # Pack ND labels along Cartesian Y to avoid overlaps.
  set nd_pairs {}
  foreach def $nd_defs {
    lassign $def name chainid colorid
    set c [_center_of_selection $mut "chain $chainid and name CA"]
    if {$c ne ""} {
      lappend nd_pairs [list $name [lindex $c $yaxis]]
    }
  }
  set nd_packed [_pack_positions $nd_pairs 22.0]
  set nd_ymap [dict create]
  foreach p $nd_packed {
    lassign $p key pos
    dict set nd_ymap $key $pos
  }

  foreach def $nd_defs {
    lassign $def name chainid colorid

    mol representation NewCartoon
    mol selection "chain $chainid and protein"
    mol color ColorID $colorid
    mol material Opaque
    mol addrep $mut

    set c [_center_of_selection $mut "chain $chainid and name CA"]
    if {$c ne ""} {
      set lpos $c
      set lpos [_set_axis_component $lpos $xaxis $label_x]
      set lpos [_set_axis_component $lpos $zaxis $nd_z]
      if {[dict exists $nd_ymap $name]} {
        set lpos [_set_axis_component $lpos $yaxis [dict get $nd_ymap $name]]
      }
      _callout $mut $c $lpos $name black 1.2
    }
  }

  # Mutation site highlight (ND6 resid 64) on mutant (no spheres/VDW).
  mol representation Bonds 0.35 12.0
  mol selection "chain m and resid 64"
  mol color ColorID 4
  mol material Opaque
  mol addrep $mut

  # --- Lipid membrane (native lipids present in model) ---
  mol representation Bonds 0.20 12.0
  mol selection "resname CDL PEE PLX"
  mol color Resname
  mol material Opaque
  mol addrep $mut

  # --- FMN and NADPH ---
  mol representation Bonds 0.25 12.0
  mol selection "resname FMN"
  mol color ColorID 6
  mol material Opaque
  mol addrep $mut

  mol representation Bonds 0.25 12.0
  mol selection "resname NDP"
  mol color ColorID 2
  mol material Opaque
  mol addrep $mut

  # NADH isn't present in these coordinates; highlight the NADH binding region near FMN.
  mol representation Bonds 0.20 12.0
  mol selection "protein and within 5 of resname FMN"
  mol color ColorID 3
  mol material Opaque
  mol addrep $mut

  # --- Quinone binding cavity (approx: residues near quinone-like ligand 8Q1) ---
  mol representation Bonds 0.25 12.0
  mol selection "resname 8Q1"
  mol color ColorID 10
  mol material Opaque
  mol addrep $mut

  mol representation Bonds 0.20 12.0
  mol selection "protein and within 5 of (resname 8Q1 and chain F)"
  mol color ColorID 10
  mol material Opaque
  mol addrep $mut

  # --- Feature callouts (ABOVE the complex) ---
  set feat_defs {
    {FMN    {resname FMN}                                       "FMN"}
    {NADPH  {resname NDP}                                       "NADPH"}
    {NADH   {protein and within 5 of resname FMN}               "NADH"}
    {Q_site {resname 8Q1 and chain F}                           "Q-site"}
    {Q_cav  {protein and within 5 of (resname 8Q1 and chain F)} "quinone binding cavity"}
  }

  set feat_pairs {}
  set feat_center [dict create]
  foreach f $feat_defs {
    lassign $f key seltext label
    set c [_center_of_selection $mut $seltext]
    if {$c ne ""} {
      dict set feat_center $key $c
      lappend feat_pairs [list $key [lindex $c $yaxis]]
    } else {
      puts "NOTE: missing feature selection: $seltext"
    }
  }

  set feat_packed [_pack_positions $feat_pairs 30.0]
  set feat_ymap [dict create]
  foreach p $feat_packed {
    lassign $p key pos
    dict set feat_ymap $key $pos
  }

  foreach f $feat_defs {
    lassign $f key seltext label
    if {![dict exists $feat_center $key]} { continue }
    set c [dict get $feat_center $key]
    set lpos $c
    set lpos [_set_axis_component $lpos $xaxis $label_x]
    set lpos [_set_axis_component $lpos $zaxis $feat_z]
    if {[dict exists $feat_ymap $key]} {
      set lpos [_set_axis_component $lpos $yaxis [dict get $feat_ymap $key]]
    }
    _callout $mut $c $lpos $label black 1.0
  }

  # --- Fe-S clusters (N1a..N6b) ---
  mol representation Bonds 0.30 12.0
  mol selection "resname SF4 FES"
  mol color Resname
  mol material Opaque
  mol addrep $mut

  # Label likely cluster identities by residue IDs in 9TI4.
  set cluster_labels {
    {N1a {resname FES and chain O and resid 301 and name FE1} black}
    {N1b {resname FES and chain M and resid 803 and name FE1} black}
    {N2  {resname SF4 and chain E and resid 302 and name FE1} black}
    {N3  {resname SF4 and chain C and resid 501 and name FE1} black}
    {N4  {resname SF4 and chain M and resid 801 and name FE1} black}
    {N5  {resname SF4 and chain M and resid 802 and name FE1} black}
    {N6a {resname SF4 and chain D and resid 301 and name FE1} black}
    {N6b {resname SF4 and chain D and resid 302 and name FE1} black}
  }
  foreach item $cluster_labels {
    lassign $item label seltext colorname
    _label_at_atom $mut $seltext $label $colorname
  }

  # Re-apply after the UI initializes (helps if .vmdrc themes set black background / axes off).
  after idle _ci_force_display_prefs
  after 200 _ci_force_display_prefs

  return
}

main
