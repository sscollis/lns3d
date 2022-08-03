plot "stab-qpnc.out" u 2:(-$4) w l
set xzeroaxis
set yrange [-1:2]
set xrange [0:35]
set mxtics 5
set mytics 5
replot
