#!/bin/bash
#
# Build SLATEC on MacOS 
#
# The default DESTination is $HOME/local/slatec
#
env FC=gfortran make 
env FC=gfortran make install
exit 0
