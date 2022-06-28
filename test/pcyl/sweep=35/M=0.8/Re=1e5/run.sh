#!/bin/bash
#set -x
#
# Make new BL mesh
#
~/git/lns3d/mesh/confpc -x3 -y3 < confpc.inp
\ln -f -s grid.dat grid.xyz
#
# Interpolate from potential solution to new BL mesh
#
~/git/lns3d/mesh/interpc -x1 -p < interpc.inp
~/git/lns3d/util/npost output.R.0 
\mv wall.dat wall.pot
\mv output.q.0 pot.q
#
# Now run lns3d to get BL mean flow
#
~/git/lns3d/src/lns3d mean.nml < mean.inp
#
# Run to convergence
#
#~/git/lns3d/src/lns3d mean.nml < res.inp
#
exit 0
