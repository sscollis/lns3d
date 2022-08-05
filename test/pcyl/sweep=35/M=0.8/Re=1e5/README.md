# Parabolic cylinder

This is $\mathsf{Ma}=0.8\, \mathsf{Pr}=1\, \theta = 35^\circ$ 
flow over a parabolic cylinder based on Collis' PhD Thesis, 
Chapter 5 where $\theta$ is the sweep angle.  We start with 
the inviscid potential flow computation using `npot`.

Note that the wall is isothermal at the theoretical 
adiabatic wall temperature for a laminar flat-plate
boundary-layer.

$$ T_w = T_\infty \left[ 1 + \frac{r(\gamma - 1)}{2} 
         \mathsf{Ma}_\infty^2 \right] $$

where $\gamma = \mathsf{c_p}/\mathsf{c_v}$ is the ratio
of specific heats, taken here to be $\gamma = 1.4$ and
$r = \sqrt{\mathsf{Pr}}$.  Note that near the leading-edge,
$r = 1$ so that, in general, there will be heat transfer 
due to this BC.  However, we use $\mathsf{Pr}=1$ so that
this BC is equivalent to an adiabatic wall.

## Mean boundary layer 

This directory is for computing the mean boundary-layer flow
given the potential flow solution at $\mathsf{Re}=10^5$.

## Running

The initial setup and run is done using:
```bash
./run.sh 
```
Then, to go to convergence (which takes a while) is
accomplished using 
```bash
./rerun.sh
```
with runs for 20,000 iterations.

## Reference

There is a reference solution in the `ref` directory stored as `mean.R.0`
along with grid and metric files.  There are also boundary-layer thickness
results from `npost` and these are generally useful although some of the 
results, like the numerical value for the `edge` may be coorupted due to the
use of a conformal mesh and not a body-fitted mesh.   Results using a body
fitted mesh are available in the `lst` directory

## Linear Stability Theory

The `lst` directory takes the reference mean solution (change `run-pro.sh` if you 
wish to use a different solution) and interpolates it to a body fitted mesh
suitable for extracting profiles and thickness statistics.   This is done
using `npost` and the resulting data can be plotted using the collestion of 
`gnuplot` scripts ending in `.com`.

### Running

Use the bash scripts in the `lst` directory, possibly changing where
the mean flow solution comes from to run various linear stability calculations. See the `README.md` therein for details. 

### Notes
  1. The thickness results can be compared to those from the reference mean 
     solution so see the impact of using a conformal versus body-fitted mesh
     (which is small for most quantities)
  2. There is a slight difference in these results relative to those in 
     Collis' PhD thesis, Ch. 5.  I have not isolated these differences yet
     but they could be due to:
       1. A better mean-flow solution -- I fixed a minor error improving the
          far-field boundary conditions that removed the dependence on
          far-field boundary condition location.
       2. These have iterated to a lower residual due to the much faster
          computers available today (my laptop is 5x faster than the Cray used
          for the PhD.
       3. There could be a difference in the wall-temperature BC -- I corrected
          what looked like a discrepency in enforcing this BC.
          
## Linearized Navier-Stokes

### Short 

The `short` directory contains an LNS run over 10 crossflow wavelengths forced using a LST eigenfunction near $s=0.8$ with $\beta=35$.  See the `run.sh` script therein. 

### LNS

The `lns` directory contains a long LNS calculation from $s=5$ to $s=25$ at $\beta=35$ containing over almost 80 crossflow streamwise wavelengths forced using a LST eigenfunction. Because of this, there is a very large range in solution values from the inflow to outflow due to the exponential growth of the crossflow instability. See the `run.sh` script and `README.md` for further details.

## Point of Contact

S. Scott Collis\
flow.physics.simulation@gmail.com
