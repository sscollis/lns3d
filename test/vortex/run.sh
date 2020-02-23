#!/bin/bash
#env GFORTRAN_CONVERT_UNIT='swap' ./lns3d < input.dat
../../util/genmesh < mesh.inp
ln -s grid.dat grid.xyz
../../util/mkvortex < vort.inp
../../src/lns3d < lns3d.inp
./mpost output.R.*
exit 0
