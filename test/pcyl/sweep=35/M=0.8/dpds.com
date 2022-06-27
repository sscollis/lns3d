set log x
set xrange [0.001:5e2]
set yrange [-0.3:0.1]
set xlabel "s"
set ylabel "dp/ds"
plot "wall.dat.0" u 2:5 w l
