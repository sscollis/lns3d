#!/bin/bash
echo "Running case $1"
~/git/lns3d/mesh/confpc -x1 < confpc.inp.$1
\ln -f -s grid.dat grid.xyz
~/git/npot/src/npot < npot.inp.0
~/git/npot/src/npot < npot.inp.1
~/git/npot/src/npot < npot.inp.2
~/git/npot/src/npot < npot.inp.3
~/git/npot/src/npot < npot.inp.3
~/git/npot/src/npot < npot.inp.2
~/git/npot/src/npot < npot.inp.1
~/git/lns3d/util/npost lns.dat
cp wall.dat wall.dat.$1
exit 0
