#!/bin/bash
if [ $# -lt 2 ];
then
  echo Usage: $0 command dirs 
fi
command=$1
shift
for dir in $@
do
  echo Checkng git $command $dir
  cd $dir 
  git $command 
  cd ..
done
exit 0
