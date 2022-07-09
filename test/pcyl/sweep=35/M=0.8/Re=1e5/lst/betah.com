#
# Plot Hartree parameter:  $\beta_h$
#
plot "betah.dat" u 1:2 w l
set xrange [0.0:30.0]
set yrange [0.0:1.2]
replot

