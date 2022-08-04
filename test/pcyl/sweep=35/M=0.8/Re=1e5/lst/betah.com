#
# Plot Hartree parameter:  $\beta_h$
#
set xrange [0.0:30.0]
set yrange [0.0:1.2]
set xlabel '$s$'
set ylabel '$\beta_h(s)$'
set title 'Evolution of Hartree parameter at $\mathsf{M}=0.8$, \
$\mathsf{Re}=10^5$,$\mathsf{Pr}=1$, $\theta=35$'
set nokey
#set key notitle invert under reverse 
#Left left spacing 2 samplen 0.7
plot "betah.dat" u 1:2 w l lw 2 title '$\beta_h$'
set terminal pict2e color texarrows font "cmr,10" size 5in,3.5in; 
set output "betah.tex"; replot
set term qt; replot
clear
exit
