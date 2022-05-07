#!/bin/bash
#
# Use default base directory if not already set
#
LNS3D_DIR="${LNS3D_DIR:=../..}"
NPOT_DIR="${NPOT_DIR:=../../../npot}"
echo LNS3D base directory = $LNS3D_DIR
$LNS3D_DIR/mesh/mse < mse.inp | tee mse.log && \
ln -f -s grid.dat grid.xyz && \
$NPOT_DIR/src/npot < npot.inp | tee npot.log && \
cp lns.dat output.R.0 && \
time $LNS3D_DIR/src/lns3d < lns3d.inp | tee lns3d.log && \
./mpost output.R.*
exit $? 
