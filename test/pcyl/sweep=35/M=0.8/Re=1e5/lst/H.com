#
# Plot boundary-layer H=theta/delta
#
set xrange [0.0:30.0]
set yrange [0.3:0.38]
set xlabel '$s$'
set ylabel '$H(s)$'
set title 'Boundary Layer shape factor for $\mathsf{M}=0.8$, $\mathsf{Re}=10^5$, $\mathsf{Pr}=1$, $\theta=35^\circ$'
#set nokey
#set key notitle invert under reverse 
#Left left spacing 2 samplen 0.7
plot "delta.dat" u 1:4 w l lw 2 title '$H(s)$'
set terminal pict2e color texarrows font "cmr,10" size 5in,3.5in; 
set output "H.tex"; replot
set term qt; replot
clear
exit
