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
\mv output.R.0 mean.R.0
\cp mean.R.0 mean.dat
$LNS3D_DIR/util/npost -p -t -Wc mean.R.0 <<EOF
1 384 32 
EOF
#
# Run stability analysis (sweep in beta)
#
$STAB_DIR/stab < stab.inp
# Output eigenvector at inflow 
#
$STAB_DIR/getevec <<EOF
evec.dat
507
0
EOF
tail -n +5 space.1 > eig.pro 
#
# make disturbance IC
#
$LNS3D_DIR/util/mkeig3d < mkeig3d.inp
$LNS3D_DIR/util/lpost3d output.R.0
#
# run lns3d
#
$LNS3D_DIR/src/lns3d < lns3d.inp
#
#
#
exit 0
