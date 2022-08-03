# Linearized Navier-Stokes run

This directory runs a LNS calculation using LNS3D that inherently includes
both surface curvature and nonparallel growth.

This run uses a pre-computed reference mean solution at \$\mathsf{Ma}=0.8\, 
\mathsf{Pr}=1\, \theta = 35^\circ$ that is stored in the repository.

The pre-computed mean solution was computed on a conformal-mesh and then is
interpolated onto a body-fitted mesh here for the LNS calculation as for 
easier comparison to linear-stability theory (see the `lst` directory).

The `run.sh` script (see below) uses `stab` to compute the parallel-flow
eigenfunction to force on the inflow for frequency $\omega=0$ and spanwise
wavenumber $\beta=35$.

The streamwise wavenumber from LST for these conditions on the infow is
\[ \alpha = -2.4927746836647E+001 - 8.2432019175349E-001 \iota \]

## Running

To run the case simply begin the run with the script:
```bash
./run.sh
```
As always, you may need to check and adjust some paths for your 
environment.

To complete the run, you will need to do:
```bash
./rerun.sh
```
which will continue the run toward convergence.

S. Scott Collis\
flow.physics.simulation@gmail.com
