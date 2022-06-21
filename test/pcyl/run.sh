#!/bin/bash
if [[ $# -lt 1 ]]; then
  echo "Usage:  run.sh case args" 
  exit 2
else
  case=$1
  shift
fi
echo "Running case $case with args $@"
set -x
~/git/lns3d/mesh/confpc -x1 $@ < confpc.inp."$case"
\ln -f -s grid.dat grid.xyz
~/git/npot/src/npot < npot.inp.0
~/git/npot/src/npot < npot.inp.1
~/git/npot/src/npot < npot.inp.2
~/git/npot/src/npot < npot.inp.3
~/git/npot/src/npot < npot.inp.3
~/git/npot/src/npot < npot.inp.2
~/git/npot/src/npot < npot.inp.1
~/git/lns3d/util/npost lns.dat
cp wall.dat wall.dat."$case"
exit 0
