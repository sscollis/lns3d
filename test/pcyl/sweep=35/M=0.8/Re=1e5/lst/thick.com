#
# Plot boundary layer thickness 
#
set xrange [0.0:30.0]
set yrange [0:0.05]
set xlabel '$s$'
set ylabel offset -1 '$\delta_1(s), \delta_2(s)$'
set title 'Boundary Layer Thickness for $\mathsf{M}=0.8$, $\mathsf{Re}=10^5$, $\mathsf{Pr}=1$, $\theta=35^\circ$'
#set nokey
#set key notitle invert under reverse 
#Left left spacing 2 samplen 0.7
set terminal pict2e color texarrows font "cmr,10" size 5in,3.5in; 
set output "thick.tex"
plot "delta.dat" u 1:2 w l lw 2 title '$\delta_1$', \
delta.dat" u 1:3 w l lw 2 title '$\delta_2$'
#set term qt; replot
clear
exit
