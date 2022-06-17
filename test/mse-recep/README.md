# Acoustic Receptivity on a MSE

This is an acoustic receptivity setup for a modified super-ellipse (MSE)
following Collis PhD thesis, Section 4.4.

There are two ways to compute the solution, time-domain and 
frequency-domain. Either way, one first needs to compute the potential 
and mean-flow solutions.

## Potential and Mean

The basic problem setup is done using
```bash
./run-mean.sh
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
     (`right=3`)
  3. These are left as exercises for others to try.

As described in the Collis thesis (Section 4.4 p. 72), an auxilary
simulation is performed to compute the wave amplitude to enforce
on the top boundary consistent with the outflow sponge.  This is done
in the subdirectory `wave` and the `amp.top` file is the wave
amplitude which is specified in the namelist input file.

### Time-Domain Solution

The time domain solution can be computed from a zero-disturbance initial 
condition `$LNS3D_DIR/util/mkdist` but I have precomputed a solution that
is fully-developed and you can start from that using
```bash
cp mean.R.0 mean.dat
\rm -f output.[Rqh].*
cp real.R.0 output.R.0
$LNS3D_DIR/src/lns3d real.nml < real.inp
./lpost.sh output.R.*
```
or use the provided script
```bash
./run-real.sh
```
then use Paraview to visualize the `q` files.  Note that `LNS3D` is 
run specifying both the regular input file `real.inp` and the auxilary 
namelist input file `real.nml`.

<p align=center>
<img src=https://github.com/sscollis/lns3d/blob/master/test/mse-recep/real-v.png>
<br>Contours of vertical velocity at t=161.9</p>

### Frequency-Domain Solution

The frequency domain solution can be computed from a zero-disturbance initial 
condition `$LNS3D_DIR/util/mkdist3d` but I have precomputed a solution that
is in the steady-state and you can start from that using
```bash
cp mean.R.0 mean.dat
\rm -f output.[Rqh].*
cp complex.R.0 output.R.0
$LNS3D_DIR/src/lns3d complex.nml < complex.inp
./lpost3d.sh output.R.*
```
or use the provided script
```base
./run-cmplx.sh
```
then use Paraview to visualize the `q` files.  Note that `LNS3D` is 
run specifying both the regular input file `complex.inp` and the auxilary 
namelist input file `complex.nml`.

Also note that the `lpost3d.sh` script outputs the imaginary component 
of the solution in the PLOT3D file by executing `lpost3d.sh -i -ij` with 
the `-i` option.  The `-ij` option makes sure that all I/O uses IJ ordering 
(not the old JI ordering used on Crays).   If you want to output the 
real-part solution or a real solution remove the `-i` option (use '-r' in
lpost3d.sh) or use the `-t` option to output the solution at a particular 
time.  To see all the options that `lpost3d.sh` supports run using 
`lpost3d.sh -h`. 

Note that the imaginary component (shown below) matches the time-domain 
solution (above) quite closely indicating that they coorespond to solutons 
at nominally similar times (really nominarlly similar phases of the solution).

<p align=center>
<img src=https://github.com/sscollis/lns3d/blob/master/test/mse-recep/cmplx-v.png>
<br>Contours of the imaginary component of vertical velocity.</p>

S. Scott Collis\
flow.physics.simulation@gmail.com
