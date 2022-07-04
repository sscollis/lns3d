#!/bin/bash
LNS3D_DIR=${LNS3D_DIR:=$HOME/git/lns3d}
$LNS3D_DIR/mesh/pc -y3 < pc.inp
\ln -f -s grid.dat grid.xyz
$LNS3D_DIR/mesh/interpc -x3 -y3 o -p -ij <<EOF
../grid.dat
../output.R.36

0.01 7 .25
0.001 7 .75

EOF
$LNS3D_DIR/util/npost -p -t -Wc output.R.0 <<EOF
1 384 1
EOF
exit 0
