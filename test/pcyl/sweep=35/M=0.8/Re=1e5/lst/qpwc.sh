#!/bin/bash
STAB_DIR=$HOME/git/stab
$STAB_DIR/stab < stab-qpwc.inp
$STAB_DIR/getax -v < getax-qpwc.inp
\mv stab.out stab-qpwc.out
exit 0