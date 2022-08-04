set title "Effect of Curvature and Nonparallel flow on growth rate"
set xlabel "s"
set ylabel "Growth Rate"
set yrange [-0.5:2]
plot "ref/stab-qpnc.out" u 2:(-$4) w l t "QPNC"
replot "ref/stab-qpwc.out" u 2:(-$4) w l t "QPWC"
replot "NPsigma.dat" u 1:11 w l t "NPWC (kinetic energy)"