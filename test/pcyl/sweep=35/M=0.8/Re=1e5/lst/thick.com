plot "delta.dat" u 1:2 w l; replot "delta.dat" u 1:3 w l
set yrange [0:0.05]
set xrange [0:30]
replot