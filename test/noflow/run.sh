#!/bin/bash
../../util/genmesh < genmesh.inp
../../util/mkini < mkini.inp
\cp ic.dat output.R.0
../../src/lns3d < lns3d.inp
../../util/npost output.R.0
../../util/npost output.R.1
../../util/npost output.R.2
../../util/npost output.R.3
../../util/npost output.R.4
../../util/npost output.R.5
\ln -sf grid.dat grid.xyz
exit 0
