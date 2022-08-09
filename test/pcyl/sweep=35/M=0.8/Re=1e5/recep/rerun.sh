#!/bin/bash
#
# setup paths
#
LNS3D_DIR=${LNS3D_DIR:=$HOME/git/lns3d}
STAB_DIR=${STAB_DIR:=$HOME/git/stab}
#
# rerun lns3d
#
$LNS3D_DIR/src/lns3d < rerun.inp
#
exit 0
