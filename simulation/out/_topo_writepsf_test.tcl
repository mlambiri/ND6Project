package require topotools
set r [catch {topo writepsf "simulation/out/_tmp.psf"} err]
puts "ret=$r err=$err"
quit
