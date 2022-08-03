#!/bin/bash
STAB_DIR=$HOME/git/stab
$STAB_DIR/stab < stab-qpnc.inp
$STAB_DIR/getax -v < getax-qpnc.inp
\mv stab.out stab-qpnc.out
exit 0
