# Simple "nonlinear" acoustic wave test

This is a "nonlinear" acoustic wave test.  The domain is from 
$x \in [0,500]$; $y \in [0:100]$ and the initial condition is a 
Guassian right-running acoustic pulse positioned at $x = 250$ 
with $\sigma = 25$ and $amp=1.0e-5$.  It makes a useful test of the 
nonlinear solver.

## LNS3D
To run execut the script
```bash
./run_lns3d.sh
```
and visualize the resulting `q`-files with Paraview.
## CNS2D
To run execut the script
```bash
./run_cns2d.sh
```
and visualize the resulting `q`-files with Paraview.

Note that `cns2d` is currently not publically distributed so that you may
now have access to this code.

S. Scott Collis\
flow.physics.simulation@gmail.com
