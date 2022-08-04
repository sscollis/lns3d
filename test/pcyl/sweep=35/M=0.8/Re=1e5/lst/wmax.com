#
# Plot crossflow velocity 
#
plot "delta.dat" u 1:5 w l t "Max crossflow"
set xrange [0.0:30.0]
set yrange [-0.125:0]
set mxtics 5
set mytics 5 
set xtics 0,5,30
set ytics -0.125,0.025
replot