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
~/git/lns3d/mesh/interpc -x1 < interpc.inp
#
# Now run lns3d to get BL mean flow
#

exit 0
