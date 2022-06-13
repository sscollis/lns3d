# 2d Wave test case

This uses the linear implicit solver from LNS3d to advanced a linear wave
on a quiesent mean flow with pure Euler.

I computed the `grid.dat`, `metric.dat`, and `mean.dat` files on an SGI so
these files must be byte swapped when using on Intel (see `run.sh` which does
so for both running and postprocessing).

The base input file is `lns3d.inp` but `imp.inp` uses the 2-step second order
time advancement.

S. Scott Collis\
flow.physics.simulation@gmail.com
