plot "stab.out" u 1:(atan2($1,$2)*180/3.14) w l
set xrange [0:200]
set yrange [90:180]
set mxtics 10
set mytics 10
replot
