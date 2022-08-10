#
# Plot bump disturbance kinetic energy 
#
set xrange [0.5:2.0]
set yrange [0:8]
set xlabel '$s$'
set ylabel offset -1 '$E_k^{1/2}$'
set title 'Evolution of $E_k^{1/2}$ for $\beta=100$, $s_w=0.7$, $\sigma_w=0.01$ $\mathsf{M}=0.8$, $\mathsf{Re}=10^5$, $\mathsf{Pr}=1$, $\theta=35^\circ$'
#set nokey
#set key notitle invert under reverse 
#Left left spacing 2 samplen 0.7
set terminal pict2e color texarrows font "cmr,10" size 5in,3.5in; 
set output "dke.tex"
plot "dke.dat" u 1:2 w l lw 2 title '$E_k^{1/2}$'
#delta.dat" u 1:3 w l lw 2 title '$\delta_2$'
#set term qt; replot
clear
exit
