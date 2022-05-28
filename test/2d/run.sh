#!/bin/bash
#
# This uses input from SGI which requires byte swapping
#
./cleanup.sh
\cp -r ic.dat output.R.0
env GFORTRAN_CONVERT_UNIT='swap' ~/git/lns3d/src/lns3d < lns3d.inp
\ln -sf grid.dat grid.xyz
env GFORTRAN_CONVERT_UNIT='swap' ../../util/npost output.R.0
env GFORTRAN_CONVERT_UNIT='swap' ../../util/npost output.R.1
env GFORTRAN_CONVERT_UNIT='swap' ../../util/npost output.R.2
env GFORTRAN_CONVERT_UNIT='swap' ../../util/npost output.R.3
exit 0
