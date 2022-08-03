#!/bin/bash
#=============================================================================
#     L i n e a r i z e d   N a v i e r - S t o k e s   D r i v e r
#=============================================================================
# Author:     S. Scott Collis
# Copyright:  S. Scott Collis(c)
# Date:       08/03/2022
#=============================================================================
#
# setup paths (your might need to change these)
#
LNS3D_DIR=${LNS3D_DIR:=$HOME/git/lns3d}
STAB_DIR=${STAB_DIR:=$HOME/git/stab}
#
# Make body-fitted mesh on parabolic cylinder
#
$LNS3D_DIR/mesh/pc -y3 < pc.inp
\ln -f -s grid.dat grid.xyz
#
# SSC: Use pre.exe to compute the metrics as they seem to be 
#      incorrect coming out of pc
#
$LNS3D_DIR/pre/src/pre.exe <<EOF
0


1
EOF
#
# replace the metrics computed by pc with those computed using pre
#
\mv metric.new metric.dat
#
# now interpolate the pre-computed reference mean solution to the new
# body-fitted mesh for the LNS calculation
#
$LNS3D_DIR/mesh/interpc -x3 -y3 -o -p -ij <<EOF
../ref/grid.dat
../ref/mean.R.0

0.01 7 .25
0.001 7 .75

EOF
#
# change names to prepare for LNS run
#
\mv output.R.0 mean.R.0
\cp mean.R.0 mean.dat
#
# visualize the mean solution
#
$LNS3D_DIR/util/npost -p -t -Wc mean.R.0 <<EOF
1 1024 128 
EOF
#
# Run stability analysis to get inflow forcing 
#
$STAB_DIR/stab < stab.inp
#
# Output eigenvector at inflow (note that we grab the eigenfunction by the
# value of the eigenvalue here)
#
$STAB_DIR/getevec -v <<EOF
evec.dat
-2.4927746836647E+001  -8.2432019175349E-001
0 0
EOF
#
# You must strip off the header for mkeig3d.  Likewise, mkeig3d generates
# inflow.dat and an initial condition projected to the body-fitted mesh
# and inflow.dat must also not have an header.
#
tail -n +5 space.1 > eig.pro 
#
# make disturbance initial condition 
#
$LNS3D_DIR/util/mkeig3d < mkeig3d.inp
#
# visualize the initial condition
#
$LNS3D_DIR/util/lpost3d output.R.0
#
# run lns3d
#
$LNS3D_DIR/src/lns3d < lns3d.inp
#
# Manually execute `rerun.sh` to compute the solution for 10,000 additional
# timesteps.
#
exit 0
