#
# Plot boundary-layer H=theta/delta
#
plot "delta.dat" u 1:4 w l
set xrange [0.0:30.0]
set yrange [0.3:0.38]
replot