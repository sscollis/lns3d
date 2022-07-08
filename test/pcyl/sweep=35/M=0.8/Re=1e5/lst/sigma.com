plot "stab.out" u 1:(-$3) w l t ""
set xrange [0:200]
set yrange [-1:4]
set mxtics 10
set mytics 10
set xzeroaxis
replot
