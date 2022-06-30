# Parabolic cylinder

This is the $M=0.8$, $Pr=1$, $\theta=35^\circ$ parabolic 
cylinder case based on Collis' PhD Thesis, Chapter 5.  
This starts with the inviscid potential flow computation 
using `npot`.

## Mean boundary layer 

This directory is for computing the mean boundary-layer flow
given the potential flow solution at Re=1e5.

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

S. Scott Collis\
flow.physics.simulation@gmail.com
