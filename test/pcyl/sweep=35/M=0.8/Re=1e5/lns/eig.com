#
# Plot eigenfunction
#
set title "Spatial Eigenfunction"
set xrange [0:0.2]
set xlabel "n"
set ylabel "Eigenfunction"
plot for [j=2:11] "eig.pro" u 1:j w l t "u(".j.")"
