#!/bin/bash
gnuplot -s H.com
gnuplot -s betah.com
gnuplot -s sigma-npwc.com 
gnuplot -s sweep.com
gnuplot -s wmax.com
gnuplot -s wave-angle.com
gnuplot -s -c "pro.com" 11
gnuplot -s thick.com
gnuplot -s alpha.com
gnuplot -s sigma-beta.com 
gnuplot -s sigma-qpnc.com 
gnuplot -s sigma-qpwc.com 
gnuplot -s alpha-npwc.com 
pdflatex pcyl-plots
exit 0
