#!/bin/bash
#
# Test of interpc interpolator (assumes that NPOT solution exists) 
#
./pc-run.sh
./interpc -o < interpc.inp 
\ln -f -s grid.dat grid.xyz
\cp int.dat output.q.0
exit $?
