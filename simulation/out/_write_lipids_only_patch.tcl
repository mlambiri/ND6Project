set psf "simulation/out/complexI_9TI4_membrane.psf"
set pdb "simulation/out/complexI_9TI4_membrane_placed.pdb"
set out_psf "simulation/out/complexI_9TI4_membrane_placed_lipidsOnly.psf"
set out_pdb "simulation/out/complexI_9TI4_membrane_placed_lipidsOnly.pdb"

set molid [mol new $psf type psf waitfor all]
mol addfile $pdb type pdb waitfor all molid $molid

set sel [atomselect $molid "resname POP"]

# Write a lipid-only patch with coordinates.
# Note: PSF written by the psf molfile plugin may omit some topology sections.
animate write psf $out_psf sel $sel
animate write pdb $out_pdb sel $sel

$sel delete
quit
