# VMD visualization: Complex I + membrane patch in one scene
#
# Usage:
#   "C:\Program Files\University of Illinois\VMD2\vmd.exe" -dispdev win -e simulation/vmd_view_complexI_9TI4_with_membrane.tcl
#
# Optional args (after -args):
#   --pdb <complex.pdb>        Complex I model (default: WT heavy; includes lipids/cofactors/Fe-S)
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

proc _ci_find_repo_root {} {
  # Try to locate the repo root robustly (some VMD builds change cwd/startup paths).
  # We define the repo root as the directory that contains:
  #   output/playwright/chatgpt_botprompts/models/ci_showhide.tcl

  set script_dir [file dirname [info script]]
  set want_rel [file join "output" "playwright" "chatgpt_botprompts" "models" "ci_showhide.tcl"]

  set candidates {}
  lappend candidates [file normalize [file join $script_dir ".."]]
  lappend candidates [file normalize $script_dir]
  lappend candidates [file normalize [pwd]]
  lappend candidates [file normalize [file join [pwd] ".."]]
  if {[info exists ::env(CI_REPO_ROOT)]} {
    lappend candidates [file normalize $::env(CI_REPO_ROOT)]
  }

  # Deduplicate
  set uniq {}
  foreach c $candidates {
    if {[lsearch -exact $uniq $c] < 0} { lappend uniq $c }
  }
  set candidates $uniq

  # 1) Direct hits
  foreach root $candidates {
    if {[file exists [file join $root $want_rel]]} {
      return $root
    }
  }

  # 2) One-level-deep search (common when cwd is e.g. C:/home/workspace but repo is C:/home/workspace/<id>/...)
  foreach parent $candidates {
    foreach sub [glob -nocomplain -types d [file join $parent "*"]] {
      if {[file exists [file join $sub $want_rel]]} {
        return [file normalize $sub]
      }
    }
  }

  return ""
}

set _ci_repo_root [_ci_find_repo_root]
if {$_ci_repo_root eq ""} {
  puts "ERROR: could not locate repo root (expected: output/playwright/chatgpt_botprompts/models/ci_showhide.tcl)."
  puts "Tip: set env var CI_REPO_ROOT to your repo root (e.g. C:/home/workspace/8196399)."
  return
}

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
  - Run ci_help for common ND/lipids/cofactors toggles.
  - Background is forced to white (reapplied after idle).}
}

proc ci_mem_help {} { _ci_mem_help }

proc ci_view {mode} {
  set m [string tolower [string trim $mode]]
  switch -exact -- $m {
    complex {
      ci_hide all
      ci_show context
      ci_show nd
      ci_show lipids
      ci_show cofactors
      ci_show fes
    }
    membrane {
      ci_hide all
      ci_show membrane
    }
    both {
      ci_hide all
      ci_show context
      ci_show nd
      ci_show lipids
      ci_show cofactors
      ci_show fes
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
    pdb [file join $_ci_repo_root "output" "playwright" "chatgpt_botprompts" "models" "complexI_9TI4_WT_heavy.pdb"] \
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
  scale by 0.92

  # Apply a consistent orientation to both molecules (optional).
  if {$do_rotate} {
    set rot {{0 0 -1 0} {-1 0 0 0} {0 1 0 0} {0 0 0 1}}
    catch {molinfo $c1 set rotate_matrix $rot}
    catch {molinfo $mem set rotate_matrix $rot}
  }

  # --- Representations ---
  # Complex: match vmd_view_complexI_9TI4_basic_003.tcl (chain-colored lines + ND tubes + optional extras).

  # Context: all protein as chain-colored lines.
  ci_add_rep_register "context" $c1 {Lines} "protein" {Chain} Opaque

  # ND subunits as thicker tubes (per-chain so colors stay distinct).
  foreach def [ci_nd_defs] {
    lassign $def name chainid colorid
    set key "nd_[string tolower $name]"
    ci_add_rep_register $key $c1 {Tube 0.60 12.0} "protein and chain $chainid" [list ColorID $colorid] Opaque
  }

  # Arm representations: membrane-near vs peripheral, using proximity to lipids (or ND neighborhood fallback).
  set membrane_dist 5.0
  set arm_include_nd 0
  lassign [ci_make_arm_selections $c1 $membrane_dist] sel_arm_mem sel_arm_per
  if {!$arm_include_nd} {
    set sel_arm_mem "$sel_arm_mem and not (chain s i j r l m)"
    set sel_arm_per "$sel_arm_per and not (chain s i j r l m)"
  }
  ci_add_rep_register "arm_membrane" $c1 {Tube 0.35 12.0} $sel_arm_mem {ColorID 2} Transparent
  ci_add_rep_register "arm_peripheral" $c1 {Tube 0.35 12.0} $sel_arm_per {ColorID 6} Transparent

  # Lipids / cofactors / Fe-S (may be empty for proteinOnly models; reps are still registered for ci_show/ci_hide).
  ci_add_rep_register "lipids" $c1 {Bonds 0.20 12.0} "resname CDL PEE PLX DGT" {Resname} Opaque
  ci_add_rep_register "cofactors" $c1 {Bonds 0.25 12.0} "resname FMN NDP 8Q1 SF4 FES" {Resname} Opaque
  ci_add_rep_register "fes" $c1 {VDW 0.80 12.0} "resname SF4 FES" {Resname} Opaque

  # Membrane: connected bilayer "strands" (no spheres).
  set repid [ci_add_rep_register "mem_strands" $mem {Bonds 0.08 12.0} "all" {ColorID 2} Opaque]
  ci_register_rep "membrane" $mem $repid

  # Default visibility.
  if {$show_mode eq "complex"} {
    ci_view complex
  } elseif {$show_mode eq "membrane"} {
    ci_view membrane
  } else {
    ci_view both
  }

  after idle _ci_force_display_prefs
  after 200 _ci_force_display_prefs

  puts "Loaded. Run: ci_mem_help  (also: ci_help)"
}

main
