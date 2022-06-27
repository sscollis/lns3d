#!/bin/bash
#
# Plots conservative variables (-c) without Plot3D normalization (-np)
#
LNS3D_DIR=${LNS3D_DIR:=$HOME/git/lns3d}
for file in $@
do
  echo Processing $file
  $LNS3D_DIR/util/npost -c -np $file >> npost.log <<EOF
35
EOF
done
exit
