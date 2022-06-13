# Acoustic Receptivity

This is an acoustic receptivity setup for a modified super-ellipse following
Collis PhD thesis, Ch. 4.

There are two ways to compute the solution, time-domain and frequency-domain.
Either way, one first needs to compute the potential and mean-flow solutions.

## Potential and Mean

The basic problem setup is done using
```bash
./run.sh
```
which makes the mesh, computes the potential solution and then copies the
`mean.R.0` file to `mean.dat`.  Note that the command to compute the mean
solution using `lns3d` is therein, but current commented out as it takes 
some signifiant computation to get to a converged meanflow solution.

## Receptivity Computations

Acoustic receptivity is computing for a planar right-running acoustic wave 
with frequency $\omega = 0.5529203070318036$ that gives a wavelength in $x$
of $\lambda_x = 125$ (recall that the speed of sound is $c=10$ and the 
freestream velocity is $U_\infty = 1$ so that the wavespeed is $a=11$ so
that $\omega = 2\pi/\lambda_x*(c+U_\infty) = 0.55292$.

The receptivity calculation can be done either in time-domain or 
frequency-domain (both are shown below).  Both approaches use an inflow
sponge that enforces the acoustic wave and an outflow sponge that 
damps all disturbances.  The top BC is characteristic forcing the 
right-running acoustic wave, the wall is no-slip and isothermal at
the mean-flow temperature.  The left BC is symmetry and the right
BC is zero disturbance (okay because the sponge has damped all
disturbances at that point).

Notes:
  1. We could have used a nonreflecting BC on the right (`right=1`)
  2. One could also try the Giles nonreflecting BC on the right
     `right=3`
  3. These are left as exercises for others to try.

The following two sections show d

### Time-Domain Solution

The time domain solution can be computed from a zero-disturbance initial 
condition `$LNS3D_DIR/util/mkdist` but I have precomputed a solution that
is fully-developed and you can start from that using
```bash
cp real.R.0 output.R.0
$LNS3D_DIR/src/lns3d real.nml < real.inp
./lpost.sh output.R.*
```
then use Paraview to visualize the `q` files.

<p align=center>
<img src=https://github.com/sscollis/lns3d/blob/master/test/receptivity/real-v.png>
<br>Contours of vertical velocity at t=161.9</p>

### Frequency-Domain Solution

The frequency domain solution can be computed from a zero-disturbance initial 
condition `$LNS3D_DIR/util/mkdist3d` but I have precomputed a solution that
is in the steady-state and you can start from that using
```bash
cp complex.R.0 output.R.0
$LNS3D_DIR/src/lns3d complex.nml < complex.inp
./lpost3d.sh output.R.*
```
then use Paraview to visualize the `q` files.  Note that the `lpost3d.sh` script
outputs the imaginary component of the solution in the PLOT3D file.  Modify
the script to output the real solution or a solution as a particular time.

Note that the imaginary component (shown below) matches the time-domain solution (above)
quite closely suggesting that they coorespond to solutons a nominally similar times
(really phases of the solution).

<p align=center>
<img src=https://github.com/sscollis/lns3d/blob/master/test/receptivity/cmplx-v.png>
<br>Contours of the imaginary component of vertical velocity.</p>

S. Scott Collis\
flow.physics.simulation@gmail.com
