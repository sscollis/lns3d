#
# Plot wave angle for sweep in beta 
#
set xrange [0:200]
set yrange [90:180]
set mxtics 10
set mytics 10
set xzeroaxis
set xlabel '$\beta$'
set ylabel '$\psi$ (deg.)'
set title 'Wave-angle from quasi-parallel LST: $\mathsf{M}=0.8$, \
$\mathsf{Re}=10^5$, $\mathsf{Pr}=1$, $\theta=35^\circ$'
set nokey
#set key notitle invert under reverse 
#Left left spacing 2 samplen 0.7
set terminal pict2e color texarrows font "cmr,10" size 5in,3.5in; 
set output "wave-angle.tex"
plot "ref/stab-beta.out" u 1:(atan2($1,$2)*180/3.14) w l lw 3
#set term qt; replot
clear
exit
