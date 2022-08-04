#
# Plot profile near max crossflow 
#
plot   "profile.".ARG1 u 1:2 w l t "rho"
replot "profile.".ARG1 u 1:3 w l t "vs"
replot "profile.".ARG1 u 1:4 w l t "vn"
replot "profile.".ARG1 u 1:5 w l t "w"
replot "profile.".ARG1 u 1:6 w l t "t"
set xrange [0.0:0.03]
set yrange [-0.5:1.5]
set mxtics 10 
set mytics 5 
set xtics 0,0.01
set ytics -0.5,0.5
replot