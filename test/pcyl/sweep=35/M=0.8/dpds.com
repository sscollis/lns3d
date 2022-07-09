#
# Gnuplot script for single dp/ds plot
#
# Note that this plots the most recent wall.dat
#
# S. Scott Collis
#
set log x
set xrange [0.001:500]
set yrange [-0.3:0.1]
set xlabel "s"
set ylabel "dp/ds"
plot "wall.dat" u 2:5 w l
