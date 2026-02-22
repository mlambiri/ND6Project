package require topotools
puts "topo commands:"
puts "  [lsort [info commands topo*]]"
puts "help for topo:" 
catch {topo help} err
puts $err
catch {topo writepsf -h} err2
puts "writepsf -h: $err2"
quit
