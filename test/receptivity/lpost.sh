#!/bin/bash
#
# Plots conservative variables (-c) without Plot3D normalization (-np)
#
LNS3D_DIR=${LNS3D_DIR:=../..}
for file in $@
do
  echo Processing $file
  $LNS3D_DIR/util/lpost3d $file >> lpost.log 
done
exit
