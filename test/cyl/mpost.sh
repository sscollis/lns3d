#!/bin/bash
LNS3D_DIR="${LNS3D_DIR:=../..}"
list=$(printf '%s\n' "$@" | sort -V)
echo $list
for file in $list 
do
  echo Processing "$file"
  $LNS3D_DIR/util/lpost3d << EOF
$file
 
EOF
done
exit
