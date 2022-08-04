#
# Plot of LST comparison for $\beta=35$ 
#
set yrange [-1:2]
set xrange [0:35]
set mxtics 5
set mytics 5
set xzeroaxis
set xlabel '$s$'
set ylabel '$\sigma$'
set title 'Quasi-parallel growth-rate comparison for:  $\beta=35$, \
$\mathsf{M}=0.8$, $\mathsf{Re}=10^5$, $\mathsf{Pr}=1$, $\theta=35^\circ$'
plot "ref/stab-qpnc.out" u 2:(-$4) w l lw 3 title 'QPNC', \
"ref/stab-qpwc.out" u 2:(-$4) w l lw 3 title 'QPWC', \
"ref/npwc-sigma.out" u 1:11 w l lw 3 title 'NPWC (KE)'
set nokey
#set key notitle invert under reverse 
#Left left spacing 2 samplen 0.7
set terminal pict2e color texarrows font "cmr,10" size 5in,3.5in; 
set output "sigma-npwc.tex"; replot
set term qt; replot
clear
exit
