# Parabolic cylinder

This is the $M=0.8$ parabolic cylinder case based on Collis'
PhD Thesis, Chapter 5.  This is the inviscid potential flow
computation using `npot` on grids of several size to show
the independence of the solution wrt to freestream BC's.

This particular case uses the incompressible solution as
both the initial condition as well as the farfield BC's.

# Running

To run the baseline case, use:
```bash
./run.sh 1
```
where `1` represents the base case input file to the mesh
generator `confpc`.  

The `run.sh` script actually runs `npot` several times, restarting
with different values of $\sigma$ in order to speed convergence.

# Plotting

The wall information is available in `wall.dat.1` and you can 
plot the pressure gradient using gnuplot with
```bash
gnuplot
plot "wall.dat.1" using 1:5 with lines
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

To plot the pressure gradient for three mesh sizes use the Python
script `dpds.py`:
```bash
./dpds.py
```
yielding the following plot

![Pressure gradient](https://github.com/sscollis/lns3d/blob/master/test/pcyl/sweep=35/M=0.8/dpds.png)

## Size of domain

There is a slight difference between `Case 0` and `Case 1` but all cases 
after `1` give nearly identical wall pressure gradients.   It would appear 
that with the improved BC's used here, a domain size of 10,000 (`Case 1`)
is adequate.

S. Scott Collis\
flow.physics.simulation@gmail.com
