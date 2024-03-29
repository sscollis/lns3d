#!/bin/bash
set -e
#==============================================================================
# LNS3d simple build script
#
# S. Scott Collis
# Copyright (c)08/03/2022
#==============================================================================
#
# Slatec
#
cd slatec/src && \
make && \
make install
cd ../..
#
# LNS3d
#
cd src && \
\ln -fs gcc.mak Makefile && \
make clean && make $@ && \
#
# LNS3d utilities
#
cd ../util && \
\ln -fs gcc.mak Makefile && \
make clean && make $@ && \
#make clean && make USE_NR=1 $@ && \
#
# LNS3d meshers
#
cd ../mesh && \
\ln -fs gcc.mak Makefile && \
make clean && make $@ && \
#make clean && make USE_NR=1 $@ && \
#
# LNS3d pre-processor
#
cd ../pre/src && \
\ln -fs gcc.mak Makefile && \
make clean && make $@
#
exit $? 
