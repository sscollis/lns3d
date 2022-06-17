#!/bin/bash
#
# Use default base directory if not already set
#
LNS3D_DIR="${LNS3D_DIR:=$HOME/git/lns3d}"
echo LNS3D base directory = $LNS3D_DIR
#
\cp -r mean.R.0 mean.dat
#
# Polish the frequency-domain solution
#
\rm -r output.[Rhq].*
\cp -f complex.R.0 output.R.0
time $LNS3D_DIR/src/lns3d complex.nml < complex.inp | tee complex.log && \
./lpost3d.sh output.R.*
#
exit $?
