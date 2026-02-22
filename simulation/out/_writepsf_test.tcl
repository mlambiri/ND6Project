set psf "simulation/out/complexI_9TI4_membrane.psf"
set pdb "simulation/out/complexI_9TI4_membrane_placed.pdb"
set molid [mol new $psf type psf waitfor all]
mol addfile $pdb type pdb waitfor all molid $molid
set sel [atomselect $molid "resname POP"]
set r1 [catch {animate write psf "simulation/out/_lipids_only_test.psf" sel $sel} err1]
puts "animate write psf ret=$r1 err=$err1"
set r2 [catch {animate write pdb "simulation/out/_lipids_only_test.pdb" sel $sel} err2]
puts "animate write pdb ret=$r2 err=$err2"
$sel delete
quit
