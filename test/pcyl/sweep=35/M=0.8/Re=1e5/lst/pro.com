#
# Plot profile near max crossflow 
#
# This script takes an argument to choose which profile to plot.
#
# You run this with:  call "pro.com" <number>
#
set xrange [0.0:0.03]
set yrange [-0.5:1.5]
set mxtics 10 
set mytics 5 
set xtics 0,0.01
set ytics -0.5,0.5
set xlabel '$n$'
#set ylabel 'Profile'
set title 'Boundary layer profile at station '.ARG1.' for $\mathsf{M}=0.8$, $\mathsf{Re}=10^5$, $\mathsf{Pr}=1$, $\theta=35^\circ$'
plot 'profile.'.ARG1 u 1:2 w l lw 3 t '$\rho$', \
'profile.'.ARG1 u 1:3 w l lw 3 t '$v_s$', \
'profile.'.ARG1 u 1:4 w l lw 3 t '$v_n$', \
'profile.'.ARG1 u 1:5 w l lw 3 t '$w$', \
'profile.'.ARG1 u 1:6 w l lw 3 t '$T$'
set terminal pict2e color texarrows font "cmr,10" size 5in,3.5in;
set output "pro.tex"; replot
set term qt; replot
clear
exit
