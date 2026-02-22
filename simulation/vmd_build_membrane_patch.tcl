# VMD helper: build a bilayer patch sized to a Complex I system
#
# Requirement: membrane size = (complex dims + 50 Å) in X/Y,
# i.e. a +25 Å margin on each side (default).
#
# Usage:
#   vmd -dispdev text -e simulation/vmd_build_membrane_patch.tcl -args --pdb system.pdb
#   vmd -dispdev text -e simulation/vmd_build_membrane_patch.tcl -args --pdb system.pdb --margin 25 --lipid POPC --out simulation/out/membrane
#
# Notes:
# - Assumes membrane plane is XY and normal is Z in the input coordinates.
# - Builds a *standalone* membrane patch (PSF+PDB). Embedding/solvation is left to your workflow.

proc _sh_parse_bool {s} {
  set t [string tolower [string trim $s]]
  if {[lsearch -exact {1 true yes on} $t] >= 0} { return 1 }
  if {[lsearch -exact {0 false no off} $t] >= 0} { return 0 }
  return -1
}

proc _usage {} {
  puts ""
  puts "vmd_build_membrane_patch.tcl"
  puts "  vmd -dispdev text -e simulation/vmd_build_membrane_patch.tcl -args --pdb system.pdb"
  puts ""
  puts "Options (after -args):"
  puts "  --pdb <file.pdb>     Input PDB (system or protein-only)"
  puts "  --psf <file.psf>     Optional PSF (if you want to load with topology)"
  puts "  --sel <selection>    Atom selection for extent calculation (default excludes water/ions/lipids)"
  puts "  --margin <A>         Margin (A) per side (default: 25)"
  puts "  --lipid <name>       Lipid name for VMD membrane plugin (default: POPC)"
  puts "  --out <prefix>       Output prefix (default: membrane_patch)"
  puts "  --place <0|1>        Translate patch to protein XY center (default: 1)"
  puts "  --z <A>              Extra Z shift applied when placing (default: 0)"
  puts "  --help               Print help and quit"
  puts ""
}

proc _parse_args {argc argv} {
  set default_sel "not (resname HOH WAT TIP TIP3 TP3 SOD CLA POT CAL MG ZN NA CL POP TYC CDL PEE PLX DGT)"
  set opts [dict create \
    help 0 \
    pdb "" \
    psf "" \
    sel $default_sel \
    margin 25.0 \
    lipid "POPC" \
    out "membrane_patch" \
    place 1 \
    z 0.0 \
  ]

  for {set i 0} {$i < $argc} {incr i} {
    set a [lindex $argv $i]
    switch -exact -- $a {
      -h -
      --help {
        dict set opts help 1
      }
      --pdb {
        incr i
        dict set opts pdb [lindex $argv $i]
      }
      --psf {
        incr i
        dict set opts psf [lindex $argv $i]
      }
      --sel {
        incr i
        dict set opts sel [lindex $argv $i]
      }
      --margin {
        incr i
        set v [lindex $argv $i]
        if {[catch {expr {double($v)}} dv]} {
          puts "WARNING: --margin expects a number (got: $v); using 25.0"
          set dv 25.0
        }
        dict set opts margin $dv
      }
      --lipid {
        incr i
        dict set opts lipid [lindex $argv $i]
      }
      --out {
        incr i
        dict set opts out [lindex $argv $i]
      }
      --place {
        incr i
        set b [_sh_parse_bool [lindex $argv $i]]
        if {$b < 0} {
          puts "WARNING: --place expects 0/1 (got: [lindex $argv $i]); using 1"
          set b 1
        }
        dict set opts place $b
      }
      --z {
        incr i
        set v [lindex $argv $i]
        if {[catch {expr {double($v)}} dv]} {
          puts "WARNING: --z expects a number (got: $v); using 0.0"
          set dv 0.0
        }
        dict set opts z $dv
      }
      default {
        if {[string match "-*" $a]} {
          puts "WARNING: unknown option: $a (ignored)"
        } elseif {[dict get $opts pdb] eq ""} {
          # allow first positional as pdb
          dict set opts pdb $a
        }
      }
    }
  }

  return $opts
}

proc main {} {
  global argc argv

  set opts [_parse_args $argc $argv]
  if {[dict get $opts help]} {
    _usage
    quit
  }

  set pdb_file [dict get $opts pdb]
  if {$pdb_file eq ""} {
    puts "ERROR: missing --pdb <file.pdb>"
    _usage
    quit
  }
  set psf_file [dict get $opts psf]
  set seltext  [dict get $opts sel]
  set margin   [dict get $opts margin]
  set lipid    [dict get $opts lipid]
  set outpref  [dict get $opts out]
  set do_place [dict get $opts place]
  set z_shift  [dict get $opts z]

  if {$psf_file ne ""} {
    set molid [mol new $psf_file type psf waitfor all]
    mol addfile $pdb_file type pdb waitfor all molid $molid
  } else {
    set molid [mol new $pdb_file type pdb waitfor all]
  }

  set sel [atomselect $molid $seltext]
  if {[$sel num] < 1} {
    puts "ERROR: selection '$seltext' is empty."
    $sel delete
    quit
  }

  set mm [measure minmax $sel]
  set minv [lindex $mm 0]
  set maxv [lindex $mm 1]
  $sel delete

  set dx [expr {[lindex $maxv 0] - [lindex $minv 0]}]
  set dy [expr {[lindex $maxv 1] - [lindex $minv 1]}]

  set patch_x [expr {$dx + 2.0*$margin}]
  set patch_y [expr {$dy + 2.0*$margin}]

  puts "Complex extents from selection '$seltext':"
  puts [format "  dx=%.3f A  dy=%.3f A" $dx $dy]
  puts [format "Membrane patch target (margin=%.3f A per side):" $margin]
  puts [format "  x=%.3f A  y=%.3f A" $patch_x $patch_y]

  if {[catch {package require membrane} err]} {
    puts "ERROR: could not load VMD membrane plugin (package 'membrane')."
    puts "  $err"
    quit
  }

  puts "Building membrane with lipid '$lipid'..."
  membrane -l $lipid -x $patch_x -y $patch_y -o $outpref

  puts "Wrote:"
  puts "  ${outpref}.psf"
  puts "  ${outpref}.pdb"

  if {$do_place} {
    set cx [expr {0.5*([lindex $minv 0] + [lindex $maxv 0])}]
    set cy [expr {0.5*([lindex $minv 1] + [lindex $maxv 1])}]
    set shift [list $cx $cy $z_shift]

    if {[file exists "${outpref}.psf"] && [file exists "${outpref}.pdb"]} {
      set memid [mol new "${outpref}.psf" type psf waitfor all]
      mol addfile "${outpref}.pdb" type pdb waitfor all molid $memid
      set memsel [atomselect $memid "all"]
      $memsel moveby $shift
      $memsel writepdb "${outpref}_placed.pdb"
      $memsel delete
      puts "Placed membrane (XY-centered on protein; +Z shift=$z_shift A):"
      puts "  ${outpref}_placed.pdb"
    } else {
      puts "WARNING: missing ${outpref}.psf/.pdb; cannot place membrane."
    }
  }

  quit
}

main
