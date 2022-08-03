# Linearized Navier-Stokes run

This directory runs a LNS calculation using LNS3D that inherently includes
both surface curvature and nonparallel growth.

This run uses a pre-computed reference mean solution at 
$\mathsf{Ma}=0.8$, $\mathsf{Pr}=1$, $\theta = 35^\circ$ 
that is stored in the repository.

The pre-computed mean solution was computed on a conformal-mesh and then is
interpolated onto a body-fitted mesh here for the LNS calculation as for 
easier comparison to linear-stability theory (see the `lst` directory).

The `run.sh` script (see below) uses `stab` to compute the parallel-flow
eigenfunction to force on the inflow for frequency $\omega=0$ and spanwise
wavenumber $\beta=35$.

The streamwise wavenumber from LST for these conditions on the infow is
$$\alpha = -24.927746836647 - 0.82432019175349 i$$

Note that the inflow conditions are now specified fully with the namelist 
input file `inflow.nm`.

The original eigenfunction, computed by `stab` can be visualized using Gnuplot
with
```bash
gnuplot
load 'eig.com'
```
and the inflow disturbance, as projected by `mkeig3d` to the global LNS 
$(x,y,z)$ coordinates using `inflow.com`.

## Running

To fully setup and begin the LNS3d calculation, execute the script:
```bash
./run.sh
```
As always, you may need to check and adjust some paths for your 
environment.

To complete the run, you will need to do:
```bash
./rerun.sh
```
which will continue the run toward convergence.  Note that this is a 
steady-state solution and that there are sponges on the outflow and top
boundaries.  Also, the `rerun.inp` input file (used by `rerun.sh`) turns 
on the fourth-order smoother to remove some minor node-to-node oscillations.

S. Scott Collis\
flow.physics.simulation@gmail.com
