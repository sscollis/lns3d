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

## Time-Domain Solution

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

## Frequency-Domain Solution

The frequency domain solution can be computed from a zero-disturbance initial 
condition `$LNS3D_DIR/util/mkdist3d` but I have precomputed a solution that
is in the steady-state and you can start from that using
```bash
cp complex.R.0 output.R.0
$LNS3D_DIR/src/lns3d complex.nml < complex.inp
./lpost3d.sh output.R.*
```
then use Paraview to visualize the `q` files.

<p align=center>
<img src=https://github.com/sscollis/lns3d/blob/master/test/receptivity/cmplx-v.png>
<br>Contours of the imaginary component of vertical velocity.</p>

S. Scott Collis\
low.physics.simulation@gmail.com
