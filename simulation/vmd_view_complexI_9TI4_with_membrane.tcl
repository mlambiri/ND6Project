# VMD visualization: Complex I + membrane patch in one scene
#
# Usage:
#   "C:\Program Files\University of Illinois\VMD2\vmd.exe" -dispdev win -e simulation/vmd_view_complexI_9TI4_with_membrane.tcl
#
# Optional args (after -args):
#   --pdb <complex.pdb>        Complex I model (default: WT proteinOnly)
#   --mem-psf <membrane.psf>   Membrane PSF (default: placed lipidsOnly)
#   --mem-pdb <membrane.pdb>   Membrane PDB (default: placed lipidsOnly)
#   --show <complex|membrane|both>  Initial view (default: both)
#   --no-recenter             Do not shift molecules near origin
#   --no-rotate               Do not apply a fixed rotate_matrix for a consistent view
#   --help
#
# Interactive (VMD Tk Console):
#   ci_mem_help
#   ci_view complex
#   ci_view membrane
#   ci_view both

set _ci_repo_root [file normalize [file join [file dirname [info script]] ".."]]
set _ci_models_dir [file join $_ci_repo_root "output" "playwright" "chatgpt_botprompts" "models"]
source [file join $_ci_models_dir "ci_showhide.tcl"]
unset _ci_models_dir

proc _ci_force_display_prefs {} {
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

proc _ci_mem_help {} {
  puts {Complex I + membrane viewer

Commands:
  ci_view complex    ;# show Complex I only
  ci_view membrane   ;# show membrane only
  ci_view both       ;# show both

Tips:
  - Run ci_list to see registered rep groups.
  - Background is forced to white (reapplied after idle).}
}

proc ci_mem_help {} { _ci_mem_help }

proc ci_view {mode} {
  set m [string tolower [string trim $mode]]
  switch -exact -- $m {
    complex {
      ci_hide all
      ci_show complex
    }
    membrane {
      ci_hide all
      ci_show membrane
    }
    both {
      ci_hide all
      ci_show complex
      ci_show membrane
    }
    default {
      puts "ci_view expects: complex | membrane | both"
    }
  }
}

proc _parse_args {argc argv} {
  global _ci_repo_root

  set opts [dict create \
    help 0 \
    pdb [file join $_ci_repo_root "output" "playwright" "chatgpt_botprompts" "models" "complexI_9TI4_WT_heavy_proteinOnly.pdb"] \
    mem_psf [file join $_ci_repo_root "simulation" "out" "complexI_9TI4_membrane_placed_lipidsOnly.psf"] \
    mem_pdb [file join $_ci_repo_root "simulation" "out" "complexI_9TI4_membrane_placed_lipidsOnly.pdb"] \
    show "both" \
    recenter 1 \
    rotate 1 \
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
      --mem-psf {
        incr i
        dict set opts mem_psf [lindex $argv $i]
      }
      --mem-pdb {
        incr i
        dict set opts mem_pdb [lindex $argv $i]
      }
      --show {
        incr i
        dict set opts show [string tolower [string trim [lindex $argv $i]]]
      }
      --no-recenter {
        dict set opts recenter 0
      }
      --no-rotate {
        dict set opts rotate 0
      }
      default {
        # allow first positional as pdb
        if {![string match "-*" $a] && [dict get $opts pdb] eq ""} {
          dict set opts pdb $a
        }
      }
    }
  }

  return $opts
}

proc main {} {
  global argc argv _ci_repo_root

  set opts [_parse_args $argc $argv]
  if {[dict get $opts help]} {
    puts "vmd_view_complexI_9TI4_with_membrane.tcl"
    puts "  vmd -dispdev win -e simulation/vmd_view_complexI_9TI4_with_membrane.tcl"
    puts "  vmd -dispdev win -e simulation/vmd_view_complexI_9TI4_with_membrane.tcl -args --pdb <complex.pdb> --show both"
    puts ""
    puts "After load, in Tk Console: ci_mem_help"
    return
  }

  set pdb_file [file normalize [dict get $opts pdb]]
  set mem_psf  [file normalize [dict get $opts mem_psf]]
  set mem_pdb  [file normalize [dict get $opts mem_pdb]]
  set show_mode [dict get $opts show]
  set do_recenter [dict get $opts recenter]
  set do_rotate   [dict get $opts rotate]

  puts "Complex I + membrane view:"
  puts "  Complex:  $pdb_file"
  puts "  Membrane: $mem_psf"
  puts "            $mem_pdb"

  if {![file exists $pdb_file]} {
    puts "ERROR: missing complex PDB: $pdb_file"
    return
  }
  if {![file exists $mem_psf]} {
    puts "ERROR: missing membrane PSF: $mem_psf"
    return
  }
  if {![file exists $mem_pdb]} {
    puts "ERROR: missing membrane PDB: $mem_pdb"
    return
  }

  _ci_force_display_prefs

  # Reset shared show/hide registry for this scene.
  ci_reset

  # Load complex.
  set c1 [mol new $pdb_file type pdb waitfor all]
  catch {graphics $c1 delete all}
  _ci_clear_all_reps $c1

  # Load membrane.
  set mem [mol new $mem_psf type psf waitfor all]
  mol addfile $mem_pdb type pdb waitfor all molid $mem
  catch {graphics $mem delete all}
  _ci_clear_all_reps $mem

  # Recenter both molecules near origin (so resetview behaves).
  if {$do_recenter} {
    set sel_center [atomselect $c1 "protein and name CA"]
    if {[$sel_center num] > 0} {
      set pcenter [measure center $sel_center]
      set shift [list \
        [expr {-1.0*[lindex $pcenter 0]}] \
        [expr {-1.0*[lindex $pcenter 1]}] \
        [expr {-1.0*[lindex $pcenter 2]}] \
      ]
      set sel_c1_all [atomselect $c1 "all"]
      $sel_c1_all moveby $shift
      $sel_c1_all delete

      set sel_mem_all [atomselect $mem "all"]
      $sel_mem_all moveby $shift
      $sel_mem_all delete
    }
    $sel_center delete
  }

  display resetview

  # Apply a consistent orientation to both molecules (optional).
  if {$do_rotate} {
    set rot {{0 0 -1 0} {-1 0 0 0} {0 1 0 0} {0 0 0 1}}
    catch {molinfo $c1 set rotate_matrix $rot}
    catch {molinfo $mem set rotate_matrix $rot}
  }

  # --- Representations ---
  # Complex: chain-colored cartoon + ND tubes.
  set repid [ci_add_rep_register "context" $c1 {NewCartoon} "protein" {Chain} Opaque]
  ci_register_rep "complex" $c1 $repid

  foreach def [ci_nd_defs] {
    lassign $def name chainid colorid
    set key "nd_[string tolower $name]"
    set repid [ci_add_rep_register $key $c1 {Tube 0.60 12.0} "protein and chain $chainid" [list ColorID $colorid] Opaque]
    ci_register_rep "complex" $c1 $repid
  }

  # Membrane: phosphorus atoms as semi-transparent spheres (fast + shows bilayer).
  set repid [ci_add_rep_register "mem_p" $mem {VDW 1.00 12.0} "name P" {ColorID 2} Transparent]
  ci_register_rep "membrane" $mem $repid

  # Default visibility.
  ci_hide all
  ci_show complex
  ci_show membrane
  if {$show_mode eq "complex"} {
    ci_view complex
  } elseif {$show_mode eq "membrane"} {
    ci_view membrane
  } else {
    ci_view both
  }

  after idle _ci_force_display_prefs
  after 200 _ci_force_display_prefs

  puts "Loaded. Run: ci_mem_help"
}

main
