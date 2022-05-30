#!/bin/bash
LNS3D_DIR="${LNS3D_DIR:=../../..}"
for file in $@
do
  echo Processing "$file"
  $LNS3D_DIR/util/lpost3d << EOF
$file
 
EOF
done
exit
