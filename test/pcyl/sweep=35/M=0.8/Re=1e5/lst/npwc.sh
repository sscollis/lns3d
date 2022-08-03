#!/bin/bash
#
# setup paths
#
LNS3D_DIR=${LNS3D_DIR:=$HOME/git/lns3d}
STAB_DIR=${STAB_DIR:=$HOME/git/stab}
SHOOT_DIR=${SHOOT_DIR:=$HOME/git/shoot}
#
# make mean LST grid and mean flow
#
$LNS3D_DIR/util/npost -l output.R.0
#
# make the metrics for LST grid
#
$LNS3D_DIR/pre/src/pre.exe <<EOF
0
lstx.dat
lstm.dat
1
EOF
#
# Compute nonparallel correction with surface curvature
#
$SHOOT_DIR/shoot-2d.exe < npwc.inp
#
exit $?
