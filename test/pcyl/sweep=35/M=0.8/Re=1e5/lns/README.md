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
$$\alpha_{LST} = -24.927746836647 - 0.82432019175349 i$$

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

## Initial Run

To fully setup and begin the LNS3d calculation, execute the script:
```bash
./run.sh
```
As always, you may need to check and adjust some paths for your 
environment.

### Visualization

The initial transient can be visualized using `Paraview` by making
local `Plot3D` files
```bash
./lpost3d.sh output.R.*
```
Then load the `grid.xyz` and `output.q.*` files into `Paraview`.

## Complete the run

To complete the run, you will need to do:
```bash
./rerun.sh
```
which will continue the run toward convergence.  Note that this is a 
steady-state solution and that there are sponges on the outflow and top
boundaries.  Also, the `rerun.inp` input file (used by `rerun.sh`) turns 
on the fourth-order smoother to remove some minor node-to-node oscillations.

From this solution, you will likely want to extract growth-rate informaton
to compare to LST.   To do so, you will run the `$LNS_DIR/util/lpost3d` utility
directly with options
```bash
$LNS3D_DIR/util/lpost3d -g -m output.R.#
```
where the particular output file that you use is converged to the steady-state
solution.   I believe that executing `rerun.sh` once is likely close enough, 
but you may want to experiment.

*Figure 1* shows a super-position of the LNS run on top of the mean-flow
for these conditions to give a sense of the two domains and different computations. 

<p align=center>
<img src=https://github.com/sscollis/lns3d/blob/master/test/pcyl/sweep=35/M=0.8/Re=1e5/lns/v-mean-cf.png>
<br><b>Figure 1:</b>Meanflow with long crossflow computation superimposed.  Contours of vertical velocity.</p>

*Figure 2* shows a closeup of the crossflow instability growing in the domain.  Note that due to exponential growth, there is a high dynamic range in these results and the contours are saturated. 

<p align=center>
<img src=https://github.com/sscollis/lns3d/blob/master/test/pcyl/sweep=35/M=0.8/Re=1e5/lns/v-cf.png>
<br><b>Figure 2:</b>Closeup of the crossflow instability for the LNS calculation on a body-fitted mesh.  Contours of vertical velocity.</p>

S. Scott Collis\
flow.physics.simulation@gmail.com
