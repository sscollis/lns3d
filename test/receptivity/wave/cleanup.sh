#!/bin/bash
\rm -f output.[Rqhd].* && \
\rm -f grid.xyz && \
\rm -f amp.dat coord.dat echo.dat fort.* *.log mag.dat sponge*.dat
\rm -f grid.dat ic.dat mean.dat metric.dat body.dat
exit $?
