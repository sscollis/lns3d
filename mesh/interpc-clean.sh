#!/bin/bash
\rm -f body.dat coord.dat wall.dat field.mean top.dat output.[Rq].0 grid.dat \
metric.dat grid.xyz *.pot int.dat \
fort.*
exit $?
