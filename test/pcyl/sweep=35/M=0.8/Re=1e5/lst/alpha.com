#
# Plot alpha for sweep in beta 
#
set xrange [0:200]
set yrange [-200:0]
set mxtics 10
set mytics 10
set xzeroaxis
set xlabel '$\beta$'
set ylabel '$\alpha$'
set title 'Quasi-parallel LST chordwise wavenumber: $\mathsf{M}=0.8$, $\mathsf{Re}=10^5$, $\mathsf{Pr}=1$, $\theta=35^\circ$'
plot "ref/stab-beta.out" u 1:($2) w l lw 3
set nokey
#set key notitle invert under reverse 
#Left left spacing 2 samplen 0.7
set terminal pict2e color texarrows font "cmr,10" size 5in,3.5in; 
set output "alpha.tex"; replot
set term qt; replot
