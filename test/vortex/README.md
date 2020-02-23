# Vortex Rebound

## Background

This is a pair of compressible vortices that rebound off a wall.
At Re=200, Ma=0.5 and Pr=1.

## How to run

The script `run.sh` shows the typical workflow and assumes that
the solver `lns3d` and utilities are located in their standard
locations relative to this directory.

To execute this example, simply type:

    ./run.sh

`mpost` runs the `lns3d` nonlinear analysis post-processor `npost` 
on all `output.R.*` files to create `output.q.*` files in Plot3d
format.

Note that Paraview prefers the Plot3d grid file to have the `xyz`
extention and the q-files to have a `q` extension.

`mkgif.sh` uses Imagemagik::convert to created an animated gif
of the `jpeg` files output from a Paraview animation.

`cleanup` is a script that returns the directory to its original
condition, ready to be run again.

## Input files

Input file     |   Description
---------------|------------------------------------------------------------
`lns3d.inp`    | Main input for LNS3d solver
`mesh.inp`     | Input for genmesh in `lns3d/utils`
`vort.inp`     | Input for `mkvortex` in `lns3d/utils' creating vortex pair

---

S. Scott Collis\
Sun Feb 23 08:02:03 MST 2020
