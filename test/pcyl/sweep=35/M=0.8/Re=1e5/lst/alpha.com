plot "stab-beta.out" u 1:($2) w l t ""
set xrange [0:200]
set yrange [-200:0]
set mxtics 10
set mytics 10
set xzeroaxis
replot
