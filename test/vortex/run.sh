#!/bin/bash
LNS3D_DIR="${LNS3D_DIR:=../..}"
echo LNS3D base directory = $LNS3D_DIR
$LNS3D_DIR/util/genmesh < mesh.inp && \
ln -f -s grid.dat grid.xyz && \
$LNS3D_DIR/util/mkvortex < vort.inp && \
$LNS3D_DIR/src/lns3d < lns3d.inp && \
./mpost output.R.*
exit 0
