#
# Plot crossflow velocity 
#
set xrange [0.0:30.0]
set yrange [-0.125:0]
set mxtics 5
set mytics 5 
set xtics 0,5,30
set ytics -0.125,0.025
set xlabel '$s$'
set ylabel offset -1 '$w_{max}(s)$'
set title 'Maximum crossflow velocity for $\mathsf{M}=0.8$, \
$\mathsf{Re}=10^5$,$\mathsf{Pr}=1$, $\theta=35$'
set nokey
plot "delta.dat" u 1:5 w l lw 2 title 'Max crossflow velocity'
set terminal pict2e color texarrows font "cmr,10" size 5in,3.5in; 
set output "wmax.tex"; replot
set term qt; replot
clear
exit
