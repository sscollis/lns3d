#!/bin/bash
set -e
#
# Execute git commands in a list of directories
#
# Note:  make sure to use single quotes 'command command' for multi word
# commands
#
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
