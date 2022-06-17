#!/bin/bash
#
# Use default base directory if not already set
#
LNS3D_DIR="${LNS3D_DIR:=$HOME/git/lns3d}"
echo LNS3D base directory = $LNS3D_DIR
#
# Run one period of the time-domain solution
#
\rm output.[Rqh].*
\cp -f real.R.0 output.R.0
time $LNS3D_DIR/src/lns3d real.nml < real.inp | tee real.log && \
./lpost.sh output.R.*
#
exit $?
