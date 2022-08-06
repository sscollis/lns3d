# Short 

This is a LNS calculation done on a relatively short domain as a simple test case for LNS3D. 

The conditions are similar to the full `lns` case so see the README therein. 

You can use the `setup.sh` script (correcting paths as needed) to setup and
`run.sh` to complete the run. 

# Plotting

Post-process the `output.R.*` files with `lpost3d.sh` and visualize the 
`grid.xyz` and `output.q.*` using Paraview.

<p align=center>
<img src=https://github.com/sscollis/lns3d/blob/master/test/pcyl/sweep=35/M=0.8/Re=1e5/short/mean-cf.png>
<br>Meanflow with short crossflow computation superimposed.  Contours of
streamwise velocity.</p>

S. Scott Collis\
flow.physics.simulation@gmail.com
