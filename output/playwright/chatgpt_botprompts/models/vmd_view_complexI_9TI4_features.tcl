# VMD visualization: single 9TI4 Complex I model with ND subunits, cofactors, lipids, and key sites marked
#
# Usage (defaults to WT full heavy model):
#   vmd -e vmd_view_complexI_9TI4_features.tcl
#
# Usage (explicit file):
#   vmd -e vmd_view_complexI_9TI4_features.tcl -args /path/to/model.pdb

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

proc _pick_membrane_normal_axis {molid} {
  # Heuristic: pick the axis with smallest extent for lipids' P atoms (bilayer thickness).
  # Fallback: ND1..ND6 CA atoms.
  set sel [atomselect $molid "resname CDL PEE PLX and (name P or element P)"]
  if {[$sel num] < 1} {
    $sel delete
    set sel [atomselect $molid "protein and (chain s i j r l m) and name CA"]
  }
  set mm [measure minmax $sel]
  $sel delete
  set minv [lindex $mm 0]
  set maxv [lindex $mm 1]
  set ext0 [expr {[lindex $maxv 0] - [lindex $minv 0]}]
  set ext1 [expr {[lindex $maxv 1] - [lindex $minv 1]}]
  set ext2 [expr {[lindex $maxv 2] - [lindex $minv 2]}]
  set axis 0
  set best $ext0
  if {$ext1 < $best} { set axis 1; set best $ext1 }
  if {$ext2 < $best} { set axis 2; set best $ext2 }
  return $axis
}

proc _set_axis_component {xyz axis value} {
  set out $xyz
  lset out $axis $value
  return $out
}

proc _label_outside {molid center label colorname offset} {
  # Place label outside the protein (reduces overlap), with a leader line.
  if {$center eq ""} { return }
  set pcenter [_center_of_selection $molid "protein and name CA"]
  if {$pcenter eq ""} { set pcenter $center }
  set dx [expr {[lindex $center 0] - [lindex $pcenter 0]}]
  set dy [expr {[lindex $center 1] - [lindex $pcenter 1]}]
  set dz [expr {[lindex $center 2] - [lindex $pcenter 2]}]
  set norm [expr {sqrt($dx*$dx + $dy*$dy + $dz*$dz)}]
  if {$norm < 1e-6} { set norm 1.0 }
  set ux [expr {$dx / $norm}]
  set uy [expr {$dy / $norm}]
  set uz [expr {$dz / $norm}]
  set lpos [list \
    [expr {[lindex $center 0] + $offset*$ux}] \
    [expr {[lindex $center 1] + $offset*$uy}] \
    [expr {[lindex $center 2] + $offset*$uz}] \
  ]
  graphics $molid color black
  graphics $molid line $center $lpos width 2
  _graphics_text $molid $lpos $label $colorname 1.4
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
  # Callout with the *label-end* segment aligned to the Cartesian Y axis (readable in
  # the "look down X" view: Z up, Y left/right).
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

  if {$argc >= 1} {
    set pdb_file [lindex $argv 0]
  } else {
    set pdb_file [file join $script_dir "complexI_9TI4_WT_heavy.pdb"]
    if {![file exists $pdb_file]} {
      set pdb_file [file join $script_dir "complexI_9TI4_WT_heavy_proteinOnly.pdb"]
    }
  }

  mol new $pdb_file type pdb waitfor all
  set molid [molinfo top]

  # Clear any lingering OpenGL drawing objects (helps if you source the script in an existing VMD session).
  catch {graphics $molid delete all}

  _ci_force_display_prefs

	  # Recenter coordinates near the origin so resetview doesn't throw the complex into a corner
	  # (PDB coords are in an absolute reference frame).
	  set pcenter [_center_of_selection $molid "protein and name CA"]
	  if {$pcenter ne ""} {
	    set shift [list \
	      [expr {-1.0*[lindex $pcenter 0]}] \
	      [expr {-1.0*[lindex $pcenter 1]}] \
	      [expr {-1.0*[lindex $pcenter 2]}] \
	    ]
	    set sel_all [atomselect $molid "all"]
	    $sel_all moveby $shift
	    $sel_all delete
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
  if {[catch {molinfo $molid set rotate_matrix $rot} err]} {
    puts "WARNING: could not set rotate_matrix: $err"
  }

	  catch {mol delrep 0 $molid}

  # Base protein.
  mol representation NewCartoon
  mol selection "protein"
  mol color Chain
  mol material Opaque
  mol addrep $molid

  # --- ND1..ND6 marked (9TI4 chain IDs from our exported PDBs) ---
  # ND labels under the complex (low Z) and other feature labels above (high Z),
  # with label positions packed along the Cartesian Y axis for readability.
  set xaxis 0
  set yaxis 1
  set zaxis 2

  # Use the full protein box for outside placement (front/below/above).
  set sel_prot [atomselect $molid "protein and name CA"]
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
  set xmargin [expr {0.08*$xspan + 6.0}]
  set zmargin [expr {0.03*$zspan + 6.0}]
  set label_x [expr {$xmin - $xmargin}]
  set nd_z   [expr {$zmin - $zmargin}]
  set feat_z [expr {$zmax + $zmargin}]

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
    set c [_center_of_selection $molid "chain $chainid and name CA"]
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
    mol addrep $molid

    set c [_center_of_selection $molid "chain $chainid and name CA"]
    if {$c ne ""} {
      set lpos $c
      set lpos [_set_axis_component $lpos $xaxis $label_x]
      set lpos [_set_axis_component $lpos $zaxis $nd_z]
      if {[dict exists $nd_ymap $name]} {
        set lpos [_set_axis_component $lpos $yaxis [dict get $nd_ymap $name]]
      }
      _callout $molid $c $lpos $name black 1.2
    }
  }

  # Lipids (if present).
  mol representation Bonds 0.20 12.0
  mol selection "resname CDL PEE PLX"
  mol color Resname
  mol material Opaque
  mol addrep $molid

  # FMN + NADPH.
  mol representation Bonds 0.25 12.0
  mol selection "resname FMN"
  mol color ColorID 6
  mol material Opaque
  mol addrep $molid

  mol representation Bonds 0.25 12.0
  mol selection "resname NDP"
  mol color ColorID 2
  mol material Opaque
  mol addrep $molid

  # NADH site near FMN (NADH not explicitly present in the model).
  mol representation Bonds 0.20 12.0
  mol selection "protein and within 5 of resname FMN"
  mol color ColorID 3
  mol material Opaque
  mol addrep $molid

  # Quinone cavity: show 8Q1 and nearby residues.
  mol representation Bonds 0.25 12.0
  mol selection "resname 8Q1"
  mol color ColorID 10
  mol material Opaque
  mol addrep $molid

  mol representation Bonds 0.20 12.0
  mol selection "protein and within 5 of (resname 8Q1 and chain F)"
  mol color ColorID 10
  mol material Opaque
  mol addrep $molid

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
    set c [_center_of_selection $molid $seltext]
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
    _callout $molid $c $lpos $label black 1.0
  }

  # Fe-S clusters.
  mol representation Bonds 0.30 12.0
  mol selection "resname SF4 FES"
  mol color Resname
  mol material Opaque
  mol addrep $molid

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
    _label_at_atom $molid $seltext $label $colorname
  }

  # Re-apply after the UI initializes (helps if .vmdrc themes set black background / axes off).
  after idle _ci_force_display_prefs
  after 200 _ci_force_display_prefs

  return
}

main
