#
# Inflow eigenfunction
#
set title "Inflow Disturbance"
set xrange [0:0.2]
set xlabel "n"
set ylabel "Eigenfunction"
plot for [j=2:11] "inflow.dat" u 1:j w l t "u(".j.")"
