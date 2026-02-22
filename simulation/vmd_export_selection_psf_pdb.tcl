# VMD helper: export a subset of atoms (selection) to a new PSF+PDB
#
# Useful for writing a "lipids-only" membrane patch (strip waters), e.g.:
#   vmd -dispdev text -e simulation/vmd_export_selection_psf_pdb.tcl -args ^
#     --psf simulation/out/complexI_9TI4_membrane.psf ^
#     --pdb simulation/out/complexI_9TI4_membrane_placed.pdb ^
#     --sel "not water" ^
#     --out simulation/out/complexI_9TI4_membrane_placed_lipidsOnly
#
# Notes:
# - When the source structure was loaded from PSF+PDB, residue names come from
#   the PSF (e.g. POPC, TIP3), while PDB output may use 3-char residue names.

proc _usage {} {
  puts ""
  puts "vmd_export_selection_psf_pdb.tcl"
  puts "  vmd -dispdev text -e simulation/vmd_export_selection_psf_pdb.tcl -args --psf in.psf --pdb in.pdb"
  puts ""
  puts "Options (after -args):"
  puts "  --psf <file.psf>     Input PSF"
  puts "  --pdb <file.pdb>     Input PDB (coordinates)"
  puts "  --sel <selection>    Atom selection (default: 'not water')"
  puts "  --out <prefix>       Output prefix (default: selection_export)"
  puts "  --help               Print help and quit"
  puts ""
}

proc _parse_args {argc argv} {
  set opts [dict create \
    help 0 \
    psf "" \
    pdb "" \
    sel "not water" \
    out "selection_export" \
  ]

  for {set i 0} {$i < $argc} {incr i} {
    set a [lindex $argv $i]
    switch -exact -- $a {
      -h -
      --help {
        dict set opts help 1
      }
      --psf {
        incr i
        dict set opts psf [lindex $argv $i]
      }
      --pdb {
        incr i
        dict set opts pdb [lindex $argv $i]
      }
      --sel {
        incr i
        dict set opts sel [lindex $argv $i]
      }
      --out {
        incr i
        dict set opts out [lindex $argv $i]
      }
      default {
        if {[string match "-*" $a]} {
          puts "WARNING: unknown option: $a (ignored)"
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

  set psf_file [dict get $opts psf]
  set pdb_file [dict get $opts pdb]
  set seltext  [dict get $opts sel]
  set outpref  [dict get $opts out]

  if {$psf_file eq "" || $pdb_file eq ""} {
    puts "ERROR: missing --psf and/or --pdb"
    _usage
    quit
  }

  set molid [mol new $psf_file type psf waitfor all]
  mol addfile $pdb_file type pdb waitfor all molid $molid

  set sel [atomselect $molid $seltext]
  set n [$sel num]
  puts "Selection: '$seltext'  atoms=$n"
  if {$n < 1} {
    puts "ERROR: selection is empty; nothing to write."
    $sel delete
    quit
  }

  animate write psf "${outpref}.psf" sel $sel
  animate write pdb "${outpref}.pdb" sel $sel

  puts "Wrote:"
  puts "  ${outpref}.psf"
  puts "  ${outpref}.pdb"

  $sel delete
  quit
}

main
