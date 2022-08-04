#
# Plot local sweep angle 
#
set xrange [0.0:30.0]
set yrange [20:100]
set mxtics 5
set mytics 10
set xtics 0,5,30
set ytics 20,20,100
set xlabel '$s$'
set ylabel '$\theta_e(s)$ (deg.)'
set title 'Evolution of local sweep angle for $\mathsf{M}=0.8$, \
$\mathsf{Re}=10^5$,$\mathsf{Pr}=1$, $\theta=35$'
set nokey
plot "edge.dat" u 1:9 w l lw 2 t "Local sweep angle"
set terminal pict2e color texarrows font "cmr,10" size 5in,3.5in; 
set output "sweep.tex"; replot
set term qt; replot
