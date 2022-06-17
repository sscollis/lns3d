#!/bin/bash
#
# Plots conservative variables (-c) without Plot3D normalization (-np)
#
LNS3D_DIR=${LNS3D_DIR:=../../..}
for file in $@
do
  echo Processing $file
  $LNS3D_DIR/util/lpost3d -ij -f $file >> lpost3d.log 
done
exit
