#
# Plot local sweep angle 
#
plot "edge.dat" u 1:9 w l t "Local sweep angle"
set xrange [0.0:30.0]
set yrange [20:100]
set mxtics 5
set mytics 10
set xtics 0,5,30
set ytics 20,20,100
replot

