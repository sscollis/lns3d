# Linearized Navier-Stokes run

This directory runs a LNS calculation using LNS3D that inherently includes
both surface curvature and nonparallel growth.

## Running

To run the case simply use the script
```bash
./run.sh
```
but as always, you may need to check and adjust some paths for your 
environment.

To complete the run, you will need to do 
```bash
./rerun.sh
```
which will continue the run to convergence.

S. Scott Collis\
flow.physics.simulation@gmail.com
