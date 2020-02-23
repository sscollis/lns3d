# Planar acoustic wave

This is a simple one-dimensional acoustic wave test using the
2d linear solver.

    ./cleanup && ./run.sh file.inp

where `file.inp` is the lns3d input file that you with to use.
If no input file is given, `rk4.inp` is used.

Input file   |  Description
-------------|------------------------------------
`beuler.inp` |  implicit 1st order backward Euler
`mid.inp`    |  implicit 2nd order midpoint rule
`rk4.inp`    |  explicit 4th order Runge-Kutta
`multi.inp`  |  implicit 2nd order 2-step backward
`mean.inp`   |  Input to `mkmean for mean.dat
`dist.inp`   |  Input for `mkdist` for disturbance initial condition
`mesh.inp`   |  `genmesh` input for grid and metrics

To run paraview on Mac, use the following:

    env LNS3D_DIR=`pwd`/../.. open /Applications/ParaView.app --args --script=`pwd`/wave.py
