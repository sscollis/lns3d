#!/bin/bash
#
# Run script for Parabolic Cylinder
#
# S. Scott Collis
#
LNS3D=$HOME/git/lns3d/src/lns3d
NPOT=$HOME/git/npot/src/npot
NPOST=$HOME/git/lns3d/util/npost
CONFPC=$HOME/git/lns3d/mesh/confpc
#set -x
if [[ $# -lt 1 ]]; then
  echo "Usage:  run.sh case args" 
  exit 2
else
  case=$1
  shift
fi
echo "Running case $case with args: $@"
$CONFPC -x1 $@ < confpc.inp."$case"
\ln -f -s grid.dat grid.xyz
$NPOT < npot.inp.0 && \
$NPOT < npot.inp.1 && \
$NPOT < npot.inp.2 && \
$NPOT < npot.inp.3 && \
$NPOT < npot.inp.3 && \
$NPOT < npot.inp.2 && \
$NPOT < npot.inp.1
$NPOST lns.dat
cp wall.dat wall.dat."$case"
cp lns.dat.q lns.q."$case"
exit 0
