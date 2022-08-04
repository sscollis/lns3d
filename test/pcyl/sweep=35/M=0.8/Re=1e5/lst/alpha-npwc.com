#
# Plot of LST wavenumber for $\beta=35$ 
#
set yrange [-200:0]
set xrange [0.1:100]
set mxtics 10 
set mytics 10 
set log x
set xzeroaxis
set xlabel '$s$'
set ylabel '$\alpha$'
set title 'Quasi-parallel wavenumber for:  $\beta=35$, \
$\mathsf{M}=0.8$, $\mathsf{Re}=10^5$, $\mathsf{Pr}=1$, $\theta=35^\circ$'
plot "ref/stab-qpnc.out" u 2:($3) w l lw 3 title 'QPNC', \
"ref/stab-qpwc.out" u 2:($3) w l lw 3 title 'QPWC', \
"ref/npwc-alpha.out" u 1:5 w l lw 3 title 'NPWC ($v_s$)'
set key inside right bottom
#set key notitle invert under reverse 
#Left left spacing 2 samplen 0.7
set terminal pict2e color texarrows font "cmr,10" size 5in,3.5in; 
set output "alpha-npwc.tex"; replot
set term qt; replot
