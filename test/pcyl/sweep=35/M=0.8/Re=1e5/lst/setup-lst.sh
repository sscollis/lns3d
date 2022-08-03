#!/bin/bash
#
# setup paths
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
\mv metric.new metric.dat
$LNS3D_DIR/mesh/interpc -x3 -y3 -o -p -ij <<EOF
../ref/grid.dat
../ref/mean.R.0

0.01 7 .25
0.001 7 .75

EOF
#
# Extract profiles
#
#$LNS3D_DIR/util/npost -p -t -Wc output.R.0 <<EOF
#1 384 1
#EOF
#
# Make an LST mesh
#
$LNS3D_DIR/util/npost -l output.R.0
#
# Make metrics for LST mesh
#
$LNS3D_DIR/pre/src/pre.exe << EOF
0
lstx.dat
lstm.dat
1
EOF
#
exit $?
