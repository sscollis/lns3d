# Auxilary planar acoustic wave
 
This run is an auxilary run of a planar acoustic wave on a uniform mean
flow with $U_\infty = 1$ and $Ma = 1$, $Pr = 1$.  Periodic BC's are used
in the $y$-direction with an inflow and outflow sponge used on the left
and right boundaries, respectfully.   The inflow sponge enforces a 
right-running acoustic wave and the outflow sponge damps all 
disturbances. To run use:
```bash
./cleanup && ./run.sh file.inp
```
The primary output of the run is the wave amplitude `amp.dat` that 
is then used as `amp.top` for the receptivity run in the parent 
directory.

S. Scott Collis\
flow.physics.simulation@gmail.com
