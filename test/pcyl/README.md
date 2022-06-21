# Parabolic cylinder

This is the $M=0.8$ parabolic cylinder case based on Collis 
PhD Thesis, Chapter 5.  This is the inviscid potential flow
computation using `npot` on grids of several size to show
the independence of the solution wrt to freestream BC's.

This particular cases uses the incompressible solution as
both the initial condition as well as the farfield BC's.

# Running

To run the baseline case, use
```bash
./run.sh 0
```
where `0` represented the base case input file to the mesh
generator `confpc`.  

The `run.sh` script actually runs `npot` several times, restarting
with different values of $\sigma$ inorder to speed convergence.

# Plotting

The wall information is available in `wall.dat.0` and you can 
plot the pressure gradient using gnuplot with
```bash
gnuplot
plot "wall.dat.0" using 1:5 with lines
```
Note that if you want to plot in the surface coordinate `s` instead
use
```bash
gnuplot 
plot "wall.dat.0" u 2:5 w l
```
or use the command file
```bash
gnuplot
load "dpds.com"
```
It turns out that the figure 5.3 in the thesis apparently is 
mislabled and actually is plotting $\partial p/\partial s$ versus
$x$, not $s$ as shown.

S. Scott Collis\
flow.physics.simulation@gmail.com
