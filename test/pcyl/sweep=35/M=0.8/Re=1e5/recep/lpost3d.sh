#!/bin/bash
LNS3D_DIR="${LNS3D_DIR:=$HOME/git/lns3d}"
list=$(printf '%s\n' "$@" | sort -V)
echo $list
for file in $list 
do
  echo Processing "$file"
  $LNS3D_DIR/util/lpost3d << EOF >> lpost3d.log
$file
 
EOF
done
exit
