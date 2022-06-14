#!/bin/bash
#
# Uses an old SGI generated mesh and metrics to have to swap Endian
#
\ln -fs grid.dat grid.xyz
export GFORTRAN_CONVERT_UNIT='swap'
\rm -f output.[Rrqh].*
cp ic_prim.dat output.res
../../../cns2d/src/cns2d cns2d_prim.inp | tee cns2d_prim.log
./mpost.sh output.r.*
exit 0
