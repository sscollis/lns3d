#!/bin/bash
LNS3D_DIR="${LNS3D_DIR:=../..}"
for file in $@
do
  echo Processing "$file"
  $LNS3D_DIR/util/npost "$file" >> npost.log 
done
exit
