# VMD visualization: human Complex I (9TI4) basic view with CLI toggles (no labels/callouts)
#
# Goal: resemble the "chain-colored lines + highlighted membrane arm" style,
# but WITHOUT any big spheres and WITHOUT any text.
#
# Usage (defaults to ND6 M64V heavy model if present, else WT):
#   vmd -e vmd_view_complexI_9TI4_basic_003.tcl
#
# Usage (explicit file):
#   vmd -e vmd_view_complexI_9TI4_basic_003.tcl -args /path/to/model.pdb
#
# Usage (options):
#   vmd -e vmd_view_complexI_9TI4_basic_003.tcl -args [pdb] --nd ND6 --arm membrane --show context,nd,arm,lipids
#   vmd -e vmd_view_complexI_9TI4_basic_003.tcl -args --show nd,arm,lipids --nd ND3 --hide cofactors,fes
#
# Options (after -args):
#   [pdb_file]                 Optional first positional arg.
#   --pdb <file.pdb>           Explicit PDB path.
#   --nd <ND1,ND2,...|all|none> ND subunits to highlight (default: all).
#   --show <default|all|csv>   Parts to show: context,nd,arm,lipids,cofactors,fes
#                              (default: context,nd,lipids,cofactors,fes).
#   --hide <csv>               Parts to hide (removes from --show).
#   --arm <membrane|peripheral|both|off>
#                              If set (and --show not provided), also enables 'arm'.
#   --arm-include-nd <0|1>     Include ND chains in the 'arm' representation (default: 0).
#   --membrane-dist <A>        Distance cutoff for membrane arm selection (default: 5.0).
#   --help                     Print this help and return.
#
# Interactive toggles (after load, in VMD Tk Console):
#   ci_list
#   ci_show lipids
#   ci_hide lipids
#   ci_focus_nd ND6
#
# Notes:
# - Lipids/cofactors require the full heavy-atom PDB; proteinOnly models won't show them.
# - No labels/callouts are drawn.

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

proc _ci_usage {} {
  puts ""
  puts "vmd_view_complexI_9TI4_basic_003.tcl"
  puts "  vmd -e vmd_view_complexI_9TI4_basic_003.tcl"
  puts "  vmd -e vmd_view_complexI_9TI4_basic_003.tcl -args /path/to/model.pdb"
  puts "  vmd -e vmd_view_complexI_9TI4_basic_003.tcl -args --nd ND6 --arm membrane --show context,nd,arm,lipids"
  puts ""
  puts "Options:"
  puts "  --pdb <file.pdb>"
  puts "  --nd <ND1,ND2,...|all|none>"
  puts "  --show <default|all|csv parts>"
  puts "  --hide <csv parts>"
  puts "  --arm <membrane|peripheral|both|off>"
  puts "  --arm-include-nd <0|1>"
  puts "  --membrane-dist <A>"
  puts "  --help"
  puts ""
  puts "Parts: context, nd, arm, lipids, cofactors, fes"
  puts ""
}

proc _ci_split_csv {s} {
  set out {}
  foreach chunk [split $s ","] {
    foreach tok [split $chunk] {
      set tok [string tolower [string trim $tok]]
      if {$tok eq ""} { continue }
      lappend out $tok
    }
  }
  return $out
}

proc _ci_unique {lst} {
  set seen [dict create]
  set out {}
  foreach x $lst {
    if {![dict exists $seen $x]} {
      dict set seen $x 1
      lappend out $x
    }
  }
  return $out
}

proc _ci_list_remove {lst remove} {
  set out {}
  foreach x $lst {
    if {[lsearch -exact $remove $x] < 0} {
      lappend out $x
    }
  }
  return $out
}

proc _ci_parse_bool {s} {
  set t [string tolower [string trim $s]]
  if {[lsearch -exact {1 true yes on} $t] >= 0} { return 1 }
  if {[lsearch -exact {0 false no off} $t] >= 0} { return 0 }
  return -1
}

proc _ci_parse_args {argc argv} {
  set opts [dict create \
    help 0 \
    pdb_file "" \
    show_explicit 0 \
    show_raw "" \
    hide_raw {} \
    nd_raw "all" \
    arm_explicit 0 \
    arm_mode "membrane" \
    arm_include_nd 0 \
    membrane_dist 5.0 \
  ]

  set positional {}

  for {set i 0} {$i < $argc} {incr i} {
    set a [lindex $argv $i]
    switch -exact -- $a {
      -h -
      --help {
        dict set opts help 1
      }
      --pdb {
        incr i
        dict set opts pdb_file [lindex $argv $i]
      }
      --show {
        incr i
        dict set opts show_explicit 1
        dict set opts show_raw [lindex $argv $i]
      }
      --hide {
        incr i
        dict lappend opts hide_raw [lindex $argv $i]
      }
      --nd -
      --nds {
        incr i
        dict set opts nd_raw [lindex $argv $i]
      }
      --arm {
        incr i
        dict set opts arm_explicit 1
        dict set opts arm_mode [string tolower [string trim [lindex $argv $i]]]
      }
      --arm-include-nd {
        incr i
        set b [_ci_parse_bool [lindex $argv $i]]
        if {$b < 0} {
          puts "WARNING: --arm-include-nd expects 0/1 (got: [lindex $argv $i]); using 0"
          set b 0
        }
        dict set opts arm_include_nd $b
      }
      --membrane-dist {
        incr i
        set v [lindex $argv $i]
        if {[catch {expr {double($v)}} dv]} {
          puts "WARNING: --membrane-dist expects a number (got: $v); using 5.0"
          set dv 5.0
        }
        dict set opts membrane_dist $dv
      }
      default {
        if {[string match "-*" $a]} {
          puts "WARNING: unknown option: $a (ignored)"
        } else {
          lappend positional $a
        }
      }
    }
  }

  if {[dict get $opts pdb_file] eq "" && [llength $positional] >= 1} {
    dict set opts pdb_file [lindex $positional 0]
  }

  return $opts
}

proc _ci_nd_defs {} {
  # name chain colorid (9TI4 chain IDs in our exported PDBs)
  # ND1=s ND2=i ND3=j ND4=r ND5=l ND6=m
  return {
    {ND1 s 1}
    {ND2 i 3}
    {ND3 j 4}
    {ND4 r 7}
    {ND5 l 0}
    {ND6 m 9}
  }
}

proc _ci_parse_nd_list {raw} {
  set t [string toupper [string trim $raw]]
  if {$t eq "" || $t eq "ALL"} { return {ND1 ND2 ND3 ND4 ND5 ND6} }
  if {$t eq "NONE" || $t eq "OFF"} { return {} }
  set out {}
  foreach tok [_ci_split_csv $t] {
    set up [string toupper $tok]
    if {[lsearch -exact {ND1 ND2 ND3 ND4 ND5 ND6} $up] >= 0} {
      lappend out $up
    }
  }
  return [_ci_unique $out]
}

proc _ci_has_any_atoms {molid seltext} {
  set sel [atomselect $molid $seltext]
  set n [$sel num]
  $sel delete
  return [expr {$n > 0}]
}

proc _ci_add_rep {molid repArgs seltext colorArgs material} {
  eval mol representation $repArgs
  mol selection $seltext
  eval mol color $colorArgs
  if {[catch {mol material $material}]} {
    mol material Opaque
  }
  mol addrep $molid
  return [expr {[molinfo $molid get numreps] - 1}]
}

proc _ci_set_reps_visible {molid repids visible} {
  foreach repid $repids {
    catch {mol showrep $molid $repid $visible}
  }
}

proc _ci_dict_lappend {dvar key value} {
  upvar 1 $dvar d
  if {![dict exists $d $key]} {
    dict set d $key {}
  }
  dict set d $key [concat [dict get $d $key] [list $value]]
}

proc _ci_make_arm_selections {molid membrane_dist} {
  set ndchains "(chain s i j r l m)"
  if {[_ci_has_any_atoms $molid "resname CDL PEE PLX"]} {
    set sel_mem "protein and within $membrane_dist of (resname CDL PEE PLX)"
    set sel_per "protein and not (within $membrane_dist of (resname CDL PEE PLX))"
  } else {
    # Fallback for protein-only models: approximate membrane-arm neighborhood by proximity to ND chains.
    set fallback_dist 10.0
    set sel_mem "protein and within $fallback_dist of (protein and $ndchains)"
    set sel_per "protein and not (within $fallback_dist of (protein and $ndchains))"
  }
  return [list $sel_mem $sel_per]
}

proc _ci_register_rep {name repid} {
  set ::ci_reps($name) $repid
}

proc ci_list {} {
  if {![info exists ::ci_molid]} {
    puts "No Complex I reps registered."
    return
  }
  puts "Registered reps (molid=$::ci_molid):"
  foreach k [lsort [array names ::ci_reps]] {
    puts "  $k -> $::ci_reps($k)"
  }
}

proc _ci_set_group {name visible} {
  if {![info exists ::ci_molid]} { return }

  set molid $::ci_molid
  set keys {}
  set up [string toupper [string trim $name]]
  set low [string tolower [string trim $name]]

  if {[info exists ::ci_reps($low)]} {
    set keys [list $low]
  } elseif {$low eq "all"} {
    set keys [array names ::ci_reps]
  } elseif {$low eq "nd"} {
    set keys [array names ::ci_reps "nd_*"]
  } elseif {$low eq "arm"} {
    set keys [array names ::ci_reps "arm_*"]
  } elseif {[lsearch -exact {context lipids cofactors fes} $low] >= 0} {
    set keys [list $low]
  } elseif {[lsearch -exact {ND1 ND2 ND3 ND4 ND5 ND6} $up] >= 0} {
    set keys [list "nd_[string tolower $up]"]
  } else {
    puts "Unknown group: $name"
    return
  }

  foreach k $keys {
    if {![info exists ::ci_reps($k)]} { continue }
    catch {mol showrep $molid $::ci_reps($k) $visible}
  }
}

proc ci_show {name} { _ci_set_group $name 1 }
proc ci_hide {name} { _ci_set_group $name 0 }

proc ci_focus_nd {ndname} {
  set up [string toupper [string trim $ndname]]
  if {[lsearch -exact {ND1 ND2 ND3 ND4 ND5 ND6} $up] < 0} {
    puts "ci_focus_nd expects ND1..ND6 (got: $ndname)"
    return
  }
  ci_hide nd
  ci_show $up
}

proc main {} {
  global argc argv

  set opts [_ci_parse_args $argc $argv]
  if {[dict get $opts help]} {
    _ci_usage
    return
  }

  set script_dir [file dirname [info script]]
  set pdb_file [dict get $opts pdb_file]
  if {$pdb_file eq ""} {
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

  set default_parts {context nd lipids cofactors fes}
  set all_parts {context nd arm lipids cofactors fes}

  if {[dict get $opts show_explicit]} {
    set raw [string tolower [string trim [dict get $opts show_raw]]]
    if {$raw eq "" || $raw eq "default"} {
      set show_parts $default_parts
    } elseif {$raw eq "all"} {
      set show_parts $all_parts
    } else {
      set show_parts [_ci_unique [_ci_split_csv $raw]]
    }
  } else {
    set show_parts $default_parts
    if {[dict get $opts arm_explicit] && [dict get $opts arm_mode] ne "off"} {
      lappend show_parts arm
      set show_parts [_ci_unique $show_parts]
    }
  }

  foreach hide_item [dict get $opts hide_raw] {
    set hide_parts [_ci_split_csv [string tolower [string trim $hide_item]]]
    set show_parts [_ci_list_remove $show_parts $hide_parts]
  }

  set arm_mode [dict get $opts arm_mode]
  if {[lsearch -exact {none off {}} $arm_mode] >= 0} {
    set arm_mode "off"
    set show_parts [_ci_list_remove $show_parts {arm}]
  }

  set nd_list [_ci_parse_nd_list [dict get $opts nd_raw]]
  set membrane_dist [dict get $opts membrane_dist]
  set arm_include_nd [dict get $opts arm_include_nd]

  puts "Complex I view (basic_003):"
  puts "  PDB: $pdb_file"
  puts "  show parts: [join $show_parts {, }]"
  if {[llength $nd_list] > 0} {
    puts "  ND highlight: [join $nd_list {, }]"
  } else {
    puts "  ND highlight: (none)"
  }
  if {[lsearch -exact $show_parts "arm"] >= 0} {
    puts "  arm: $arm_mode (membrane_dist=$membrane_dist A, include_nd=$arm_include_nd)"
  }

  set molid [mol new $pdb_file type pdb waitfor all]
  set ::ci_molid $molid
  catch {unset ::ci_reps}
  array set ::ci_reps {}

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

  # --- Create representations (register all; visibility is set afterwards) ---

  # Context: all protein as chain-colored lines.
  set rep_context [_ci_add_rep $molid {Lines} "protein" {Chain} Opaque]
  _ci_register_rep "context" $rep_context

  # ND subunits as thicker tubes (per-chain so colors stay distinct).
  foreach def [_ci_nd_defs] {
    lassign $def name chainid colorid
    set rep [_ci_add_rep $molid {Tube 0.60 12.0} "protein and chain $chainid" [list ColorID $colorid] Opaque]
    _ci_register_rep "nd_[string tolower $name]" $rep
  }

  # Arm representations: membrane-near vs peripheral, using proximity to lipids (or ND neighborhood fallback).
  lassign [_ci_make_arm_selections $molid $membrane_dist] sel_arm_mem sel_arm_per
  if {!$arm_include_nd} {
    set sel_arm_mem "$sel_arm_mem and not (chain s i j r l m)"
    set sel_arm_per "$sel_arm_per and not (chain s i j r l m)"
  }
  set rep_arm_mem [_ci_add_rep $molid {Tube 0.35 12.0} $sel_arm_mem {ColorID 2} Transparent]
  _ci_register_rep "arm_membrane" $rep_arm_mem
  set rep_arm_per [_ci_add_rep $molid {Tube 0.35 12.0} $sel_arm_per {ColorID 6} Transparent]
  _ci_register_rep "arm_peripheral" $rep_arm_per

  # Lipids.
  set rep_lipids [_ci_add_rep $molid {Bonds 0.20 12.0} "resname CDL PEE PLX" {Resname} Opaque]
  _ci_register_rep "lipids" $rep_lipids

  # Cofactors/clusters as bonds (small, no spheres).
  set rep_cofactors [_ci_add_rep $molid {Bonds 0.25 12.0} "resname FMN NDP 8Q1 SF4 FES" {Resname} Opaque]
  _ci_register_rep "cofactors" $rep_cofactors

  # Fe-S clusters as spheres (small selection; avoids the "big spheres" issue).
  set rep_fes [_ci_add_rep $molid {VDW 0.80 12.0} "resname SF4 FES" {Resname} Opaque]
  _ci_register_rep "fes" $rep_fes

  # --- Apply initial visibility based on CLI toggles ---

  # Start with everything hidden, then enable requested parts.
  foreach k [array names ::ci_reps] {
    catch {mol showrep $molid $::ci_reps($k) 0}
  }

  if {[lsearch -exact $show_parts "context"] >= 0} { ci_show context }
  if {[lsearch -exact $show_parts "lipids"] >= 0} { ci_show lipids }
  if {[lsearch -exact $show_parts "cofactors"] >= 0} { ci_show cofactors }
  if {[lsearch -exact $show_parts "fes"] >= 0} { ci_show fes }

  if {[lsearch -exact $show_parts "nd"] >= 0} {
    ci_hide nd
    foreach nd $nd_list { ci_show $nd }
  }

  if {[lsearch -exact $show_parts "arm"] >= 0} {
    ci_hide arm
    switch -exact -- $arm_mode {
      membrane { ci_show arm_membrane }
      peripheral { ci_show arm_peripheral }
      both {
        ci_show arm_membrane
        ci_show arm_peripheral
      }
      default { }
    }
  }

  # Re-apply prefs after startup scripts/themes.
  after idle _ci_force_display_prefs_basic
  after 200 _ci_force_display_prefs_basic

  return
}

main
