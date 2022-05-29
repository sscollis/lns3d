#!/bin/bash
#
# Default input is explicit
#
INP="${1:-rk4.inp}"
LNS3D_DIR="${LNS3D_DIR:=../..}"
echo LNS3D base directory = $LNS3D_DIR
echo "  Using input file $INP"
export PATH=$LNS3D_DIR/util:$PATH
$LNS3D_DIR/util/genmesh < mesh.inp && \
ln -f -s grid.dat grid.xyz && \
$LNS3D_DIR/util/mkmean < mean.inp && \
$LNS3D_DIR/util/mkdist < dist.inp && \
$LNS3D_DIR/src/lns3d < $INP | tee lns3d.log && \
./mpost output.R.*
exit
