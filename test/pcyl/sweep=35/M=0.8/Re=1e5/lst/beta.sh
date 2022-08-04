#!/bin/bash
#
# setup paths
#
LNS3D_DIR=${LNS3D_DIR:=$HOME/git/lns3d}
STAB_DIR=${STAB_DIR:=$HOME/git/stab}
#
# Run stability analysis (sweep in beta)
#
$STAB_DIR/stab < stab-beta.inp
#
# Extract most unstable mode
#
$STAB_DIR/getab -v < getab-beta.inp
\mv stab.out stab-beta.out
#
# Output eigenvector at max instability
#
$STAB_DIR/getevec -v <<EOF
eig.14
-1.0420882789557E+002  -3.9004081609224E+000
0 0
EOF
\mv space.1 space-beta.14
exit $?