# Receptivity to surface roughness

The conditions for bump this bump receptivity computation are:
$\mathsf{M}=0.8$, $\mathsf{Re}=10^5$, $\mathsf{Pr}=1$, $\beta=100$
with bump parameters: $\sigma_w=0.01$, bump location $s_w = 0.7$ and domain 
$s\in[0.32, 1.80]$ which is equivalent to $x\in[0.05,1.0].

## Running

This case is run using first a setup of the mesh, meanflow, etc. and then
a run script.  The runscript can be executed multiple times to do restarts
but is setup to be sufficient for one run to achieve a suitable steady-state.

```bash
./setup.sh
./run.sh
```

## Results

The interplated mean solution:

![Mean Streamwise Velocity](https://github.com/sscollis/lns3d/blob/master/test/pcyl/sweep=35/M=0.8/Re=1e5/recep/u-mean.png)

The disturbance solution for the forced surface roughness:

![Disturbance Spanwise Velocity](https://github.com/sscollis/lns3d/blob/master/test/pcyl/sweep=35/M=0.8/Re=1e5/recep/v-bump.png)

The disturbance kinetic energy to be compared with Fig. 5.50 of Collis PhD
thesis on p.158 (see `docs` for a copy).

![Bump Receptivity](https://github.com/sscollis/lns3d/blob/master/test/pcyl/sweep=35/M=0.8/Re=1e5/recep/dke.png)

S. Scott Collis\
flow.physics.simulation@gmail.com
