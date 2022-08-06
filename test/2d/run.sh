#!/bin/bash
#
# This uses data pre-computed on SGI which requires byte swapping for Intel
#
INP="${1:-lns3d.inp}"
LNS3D_DIR="${LNS3D_DIR:=../..}"
echo LNS3D base directory = $LNS3D_DIR
echo "  Using input file $INP"
#export PATH=$LNS3D_DIR/util:$PATH
echo $PATH
./cleanup.sh
\cp -r ic.dat output.R.0
env GFORTRAN_CONVERT_UNIT='swap' $LNS3D_DIR/src/lns3d < lns3d.inp
\ln -sf grid.dat grid.xyz
env GFORTRAN_CONVERT_UNIT='swap' $LNS3D_DIR/util/npost output.R.0
env GFORTRAN_CONVERT_UNIT='swap' $LNS3D_DIR/util/npost output.R.1
env GFORTRAN_CONVERT_UNIT='swap' $LNS3D_DIR/util/npost output.R.2
env GFORTRAN_CONVERT_UNIT='swap' $LNS3D_DIR/util/npost output.R.3
exit 0
