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
$LNS3D_DIR/util/npost -p -t -Wc output.R.0 <<EOF
1 384 1
EOF
#
# Run stability analysis (sweep in beta)
#
$STAB_DIR/stab < stab.inp
#
# Extract most unstable mode
#
$STAB_DIR/getab < getab.inp
\rm stab.out stab-beta.out
#
# Output eigenvector at max instability
#
$STAB_DIR/getevec <<EOF
eig.14
505
0
EOF
\mv space.1 space-beta.14
exit 0
