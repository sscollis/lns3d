#!/bin/bash
if [[ $# -ne 1 ]]; then
  echo Usage:  $(basename "$0") Re
  exit 2
fi
awk -v Re=$1 'NR>2 {print sqrt(2*$1*Re), -$3*(2*$1+1)/sqrt(2*$1*Re)}' wall.dat 
exit 0
