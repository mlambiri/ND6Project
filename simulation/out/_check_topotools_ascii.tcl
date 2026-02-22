if {[catch {package require topotools} err]} {
  puts "NO topotools: $err"
} else {
  puts "topotools OK"
}
quit
