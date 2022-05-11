#!/bin/bash
ln -f -s grid.dat grid.xyz
../../util/mkdist3d < ic.inp
#env GFORTRAN_CONVERT_UNIT='swap' \
../../src/lns3d < rk4.inp
./mpost.sh output.R.*
