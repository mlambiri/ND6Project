package require membrane
puts "membrane package loaded"
catch {membrane} err
puts "membrane call result: $err"
quit
