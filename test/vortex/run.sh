#!/bin/bash
#
# Use default base directory if not 
# already set
#
LNS3D_DIR="${LNS3D_DIR:=../..}"
echo LNS3D base directory = $LNS3D_DIR
$LNS3D_DIR/util/genmesh < mesh.inp && \
ln -f -s grid.dat grid.xyz && \
$LNS3D_DIR/util/mkvortex < vort.inp && \
time $LNS3D_DIR/src/lns3d < lns3d.inp | tee lns3d.log && \
./mpost output.R.*
exit 0
