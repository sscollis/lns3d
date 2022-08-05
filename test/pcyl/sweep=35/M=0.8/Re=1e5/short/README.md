# Short 

This is a LNS calculation done on a relatively short domain as a simple test case for LNS3D. 

The conditions are similar to the full `lns` case so see the README therein. 

You can use the `run.sh` script (correcting paths as needed) to setup this run. 

There are two input files
```bash
lns3d < lns3d.inp
```
and 
```bash
lns3d < res.inp
```
for initial and final run. 

# Plotting

Post-process the `output.R.*` files with `lpost3d.sh` and visualize the `grid.xyz` and `output.a.*` using Paraview. 

S. Scott Collis\
flow.physics.simulation@gmail.com