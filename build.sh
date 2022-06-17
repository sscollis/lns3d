#!/bin/bash
cd src
\ln -fs gcc.mak Makefile
make clean && make
cd ../util
\ln -fs gcc.mak Makefile
make clean && make USE_NR=1
cd ../mesh
\ln -fs gcc.mak Makefile
make clean && make USE_NR=1
exit 0
