#!/bin/bash
#
# Uses an old SGI generated mesh and metrics to have to swap Endian
#
\ln -fs grid.dat grid.xyz
export GFORTRAN_CONVERT_UNIT='swap'
../../util/mkini << EOF
1000 1 1
0
g
1e-5 250 25
EOF
\mv output.R.0 output.res
../../../cns2d/src/cns2d cns2d_cons.inp | tee cns2d_cons.log
#./mpost.sh output.R.*
exit 0
