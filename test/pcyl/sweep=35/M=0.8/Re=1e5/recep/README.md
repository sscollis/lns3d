# Receptivity to surface roughness

The conditions for bump this bump receptivity computation are:
$\mathsf{M}=0.8$, $\mathsf{Re}=10^5$, $\mathsf{Pr}=1$, $\beta=100$
with bump parameters: $\sigma_w=0.01$, bump location $s_w = 0.7$ and domain 
$s\in[0.32, 1.80]$ which is equivalent to $x\in[0.05,1.0]$.

Notes:
  1. This corresponds to `case 6` in Table 5.2 from Collis PhD       
     (not `case 3` which was an error in the thesis)
  2. There is an outflow sponge over the last 20% of the domain, enforced
     in computational space.
  3. We use implict time advancedment with CFL=10 and first-order Euler 
     (since time accuracy is not needed).
  4. The `wall=3` boundary condition corresponds to a forced linearized
     bump (see Appendix F of Collis' PhD) and we use the version for an 
     iso-thermal wall, `wallt=1`

## Running

This case is run using first a setup of the mesh, meanflow, etc. and then
a run script.  The runscript can be executed multiple times to do restarts
but is setup to be sufficient for one run to achieve a suitable steady-state.

```bash
./setup.sh
./run.sh
```
Note that there is no check to make sure that `setup.sh` has been run, so
make sure to do so be executing `run.sh`.  This also assumes that LNS3D is
installed at `$HOME/git/lns3d`.  If your installation is elsewhere, you can
use
```bash
env LNS3D_DIR=$HOME/git/lns3d ./setup.sh
env LNS3D_DIR=$HOME/git/lns3d ./run.sh
```
to override the default value of `LNS3D_DIR` to whereever you have installed
LNS3D.  You can return the setup to its original state using `cleanup.sh`.

## Results

The interplated mean solution:

![Mean Streamwise Velocity](https://github.com/sscollis/lns3d/blob/master/test/pcyl/sweep=35/M=0.8/Re=1e5/recep/u-mean.png)

The disturbance solution for the forced surface roughness:

![Disturbance Spanwise Velocity](https://github.com/sscollis/lns3d/blob/master/test/pcyl/sweep=35/M=0.8/Re=1e5/recep/w-bump.png)

The disturbance kinetic energy to be compared with Fig. 5.50 of Collis PhD
thesis on p.158 (see `docs` for a copy).

![Bump Receptivity](https://github.com/sscollis/lns3d/blob/master/test/pcyl/sweep=35/M=0.8/Re=1e5/recep/dke.png)

## Point of contact

S. Scott Collis\
flow.physics.simulation@gmail.com
