#!/bin/bash
STAB_DIR=$HOME/git/stab
$STAB_DIR/stab < stab-s.inp
$STAB_DIR/getax -v < getax.inp
\mv stab.out stab-s.out
exit 0
