#!/bin/bash
../../mesh/cyl < cyl.inp
../../util/mkmean < mkmean.inp
../../util/mkdist3d < mkdist3d.inp
\cp output.R.0 ic.dat
\rm output.[Rqh].*
\cp ic.dat output.R.0
../../src/lns3d < lns3d.inp
./mpost.sh output.R.*
exit 0
