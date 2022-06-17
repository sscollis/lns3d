#!/bin/bash
#
# Test of pc mesh generator
#
./pc -x4 < pc.inp 
\ln -f -s grid.dat grid.xyz
exit $?
