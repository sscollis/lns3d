#!/bin/bash
#=============================================================================
#     L i n e a r i z e d   N a v i e r - S t o k e s   D r i v e r
#=============================================================================
# Author:     S. Scott Collis
# Copyright:  S. Scott Collis(c)
# Date:       08/09/2022
#
# Purpose:    Run a bump receptivity calculation similar to case in Collis'
#             PhD thesis, case 3 of Table 5.2
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
$LNS3D_DIR/util/npost -t -Wc mean.R.0
#
\rm output.[Rhq].*
$LNS3D_DIR/util/mkdist3d << EOF
Z
EOF
#
# visualize the disturbance initial condition
#
$LNS3D_DIR/util/lpost3d output.R.0
#
exit $?
