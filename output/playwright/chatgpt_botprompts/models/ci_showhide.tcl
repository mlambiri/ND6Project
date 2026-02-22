# Shared VMD Tcl helpers: register reps + toggle visibility (used by *_003 and other scripts)
#
# This file defines interactive Tk Console commands:
#   ci_list
#   ci_show <group>
#   ci_hide <group>
#   ci_focus_nd ND6
#
# And helper procs for scripts:
#   ci_reset
#   ci_register_rep <key> <molid> <repid>
#   ci_add_rep / ci_add_rep_register
#   ci_parse_cli / ci_resolve_visibility / ci_apply_visibility
#
# Groups expected by vmd_view_complexI_9TI4_basic_003.tcl:
#   context, nd (ND1..ND6), arm (arm_membrane/arm_peripheral), lipids, cofactors, fes

if {[info exists ::ci_showhide_loaded] && $::ci_showhide_loaded} {
  return
}
set ::ci_showhide_loaded 1

proc _ci_sh_split_csv {s} {
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

proc _ci_sh_unique {lst} {
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

proc _ci_sh_list_remove {lst remove} {
  set out {}
  foreach x $lst {
    if {[lsearch -exact $remove $x] < 0} {
      lappend out $x
    }
  }
  return $out
}

proc _ci_sh_parse_bool {s} {
  set t [string tolower [string trim $s]]
  if {[lsearch -exact {1 true yes on} $t] >= 0} { return 1 }
  if {[lsearch -exact {0 false no off} $t] >= 0} { return 0 }
  return -1
}

proc ci_reset {} {
  catch {unset ::ci_reps}
  array set ::ci_reps {}
}

proc ci_register_rep {name molid repid} {
  set key [string tolower [string trim $name]]
  if {$key eq ""} { return }
  if {![info exists ::ci_reps($key)]} {
    set ::ci_reps($key) {}
  }
  lappend ::ci_reps($key) [list $molid $repid]
}

proc ci_add_rep {molid repArgs seltext colorArgs material} {
  eval mol representation $repArgs
  mol selection $seltext
  eval mol color $colorArgs
  if {[catch {mol material $material}]} {
    mol material Opaque
  }
  mol addrep $molid
  return [expr {[molinfo $molid get numreps] - 1}]
}

proc ci_add_rep_register {key molid repArgs seltext colorArgs material} {
  set repid [ci_add_rep $molid $repArgs $seltext $colorArgs $material]
  ci_register_rep $key $molid $repid
  return $repid
}

proc ci_nd_defs {} {
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

proc ci_parse_nd_list {raw} {
  set t [string toupper [string trim $raw]]
  if {$t eq "" || $t eq "ALL"} { return {ND1 ND2 ND3 ND4 ND5 ND6} }
  if {$t eq "NONE" || $t eq "OFF"} { return {} }
  set out {}
  foreach tok [_ci_sh_split_csv $t] {
    set up [string toupper $tok]
    if {[lsearch -exact {ND1 ND2 ND3 ND4 ND5 ND6} $up] >= 0} {
      lappend out $up
    }
  }
  return [_ci_sh_unique $out]
}

proc _ci_sh_has_any_atoms {molid seltext} {
  set sel [atomselect $molid $seltext]
  set n [$sel num]
  $sel delete
  return [expr {$n > 0}]
}

proc ci_make_arm_selections {molid membrane_dist} {
  set ndchains "(chain s i j r l m)"
  if {[_ci_sh_has_any_atoms $molid "resname CDL PEE PLX"]} {
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

proc ci_list {} {
  if {![info exists ::ci_reps]} {
    puts "No reps registered."
    return
  }
  puts "Registered rep groups:"
  foreach k [lsort [array names ::ci_reps]] {
    puts "  $k -> $::ci_reps($k)"
  }
}

proc _ci_sh_set_group {name visible} {
  if {![info exists ::ci_reps]} { return }
  set low [string tolower [string trim $name]]
  set up  [string toupper [string trim $name]]

  set keys {}
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
    foreach pair $::ci_reps($k) {
      set molid [lindex $pair 0]
      set repid [lindex $pair 1]
      catch {mol showrep $molid $repid $visible}
    }
  }
}

proc ci_show {name} { _ci_sh_set_group $name 1 }
proc ci_hide {name} { _ci_sh_set_group $name 0 }

proc ci_focus_nd {ndname} {
  set up [string toupper [string trim $ndname]]
  if {[lsearch -exact {ND1 ND2 ND3 ND4 ND5 ND6} $up] < 0} {
    puts "ci_focus_nd expects ND1..ND6 (got: $ndname)"
    return
  }
  ci_hide nd
  ci_show $up
}

proc ci_parse_cli {argc argv} {
  set opts [dict create \
    help 0 \
    show_explicit 0 \
    show_raw "" \
    hide_raw {} \
    nd_explicit 0 \
    nd_raw "all" \
    arm_explicit 0 \
    arm_mode "membrane" \
    arm_include_nd 0 \
    membrane_dist 5.0 \
    pdb_file "" \
    wt_file "" \
    mut_file "" \
    positionals {} \
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
      --wt {
        incr i
        dict set opts wt_file [lindex $argv $i]
      }
      --mut {
        incr i
        dict set opts mut_file [lindex $argv $i]
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
        dict set opts nd_explicit 1
        dict set opts nd_raw [lindex $argv $i]
      }
      --arm {
        incr i
        dict set opts arm_explicit 1
        dict set opts arm_mode [string tolower [string trim [lindex $argv $i]]]
      }
      --arm-include-nd {
        incr i
        set b [_ci_sh_parse_bool [lindex $argv $i]]
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

  dict set opts positionals $positional
  return $opts
}

proc ci_resolve_visibility {opts {default_parts {context nd lipids cofactors fes}} {all_parts {context nd arm lipids cofactors fes}}} {

  if {[dict get $opts show_explicit]} {
    set raw [string tolower [string trim [dict get $opts show_raw]]]
    if {$raw eq "" || $raw eq "default"} {
      set show_parts $default_parts
    } elseif {$raw eq "all"} {
      set show_parts $all_parts
    } else {
      set show_parts [_ci_sh_unique [_ci_sh_split_csv $raw]]
    }
  } else {
    set show_parts $default_parts
    if {[dict get $opts arm_explicit] && [dict get $opts arm_mode] ne "off"} {
      lappend show_parts arm
      set show_parts [_ci_sh_unique $show_parts]
    }
  }

  foreach hide_item [dict get $opts hide_raw] {
    set hide_parts [_ci_sh_split_csv [string tolower [string trim $hide_item]]]
    set show_parts [_ci_sh_list_remove $show_parts $hide_parts]
  }

  set arm_mode [dict get $opts arm_mode]
  if {[lsearch -exact {none off {}} $arm_mode] >= 0} {
    set arm_mode "off"
    set show_parts [_ci_sh_list_remove $show_parts {arm}]
  }

  set nd_list [ci_parse_nd_list [dict get $opts nd_raw]]

  return [dict create \
    show_parts $show_parts \
    nd_list $nd_list \
    arm_mode $arm_mode \
    arm_include_nd [dict get $opts arm_include_nd] \
    membrane_dist [dict get $opts membrane_dist] \
  ]
}

proc ci_hide_registered {} {
  if {![info exists ::ci_reps]} { return }
  set seen [dict create]
  foreach k [array names ::ci_reps] {
    foreach pair $::ci_reps($k) {
      set molid [lindex $pair 0]
      set repid [lindex $pair 1]
      set id "${molid}:${repid}"
      if {[dict exists $seen $id]} { continue }
      dict set seen $id 1
      catch {mol showrep $molid $repid 0}
    }
  }
}

proc ci_apply_visibility {show_parts nd_list arm_mode} {
  ci_hide_registered

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
}
