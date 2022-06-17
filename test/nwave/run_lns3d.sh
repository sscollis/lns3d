#!/bin/bash
#
# Uses an old SGI generated mesh and metrics to have to swap Endian
#
\ln -fs grid.dat grid.xyz
export GFORTRAN_CONVERT_UNIT='swap'
\rm output.[Rrqh].*
cp ic_lns3d.dat output.R.0
../../src/lns3d < lns3d.inp | tee lns3d.log
./mpost.sh output.R.*
exit 0
