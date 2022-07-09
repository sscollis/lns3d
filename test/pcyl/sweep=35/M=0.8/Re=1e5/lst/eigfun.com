#
# Plot eigenvector at max instability 
#
plot   "space-beta.14" u 1:2 w l t "rho_r"
replot "space-beta.14" u 1:3 w l t "rho_i"
replot "space-beta.14" u 1:4 w l t "u_r"
replot "space-beta.14" u 1:5 w l t "u_i"
replot "space-beta.14" u 1:6 w l t "v_r"
replot "space-beta.14" u 1:7 w l t "v_i"
replot "space-beta.14" u 1:8 w l t "w_r"
replot "space-beta.14" u 1:9 w l t "w_i"
replot "space-beta.14" u 1:10 w l t "t_r"
replot "space-beta.14" u 1:11 w l t "t_i"
set xrange [0.0:0.03]
set yrange [-0.5:1.0]
set mxtics 10 
set mytics 5 
set xtics 0,0.01
set ytics -0.5,0.5
replot
