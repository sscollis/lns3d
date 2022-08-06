#!/bin/bash
#
# setup paths
#
LNS3D_DIR=${LNS3D_DIR:=$HOME/git/lns3d}
STAB_DIR=${STAB_DIR:=$HOME/git/stab}
#
# Run the case
#
$LNS3D_DIR/src/lns3d < lns3d.inp
#
exit 0
