#!/bin/bash
#
# Plots conservative variables (-c) without Plot3D normalization (-np)
#
LNS3D_DIR=${LNS3D_DIR:=$HOME/git/lns3d}
for file in $@
do
  echo Processing file $file
  $LNS3D_DIR/util/lpost3d -i -ij $file >> lpost3d.log 
done
exit
