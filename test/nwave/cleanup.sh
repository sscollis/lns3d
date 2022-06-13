#!/bin/bash
\rm -f grid.xyz
\rm output.*
\rm *.plt
\cp ic_prim.dat output.R.0
\cp ic_prim.dat output.r.000000
\cp ic_prim.dat output.res
\rm *.log echo.dat wall.dat
exit 0
