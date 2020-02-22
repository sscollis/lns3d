# LNS3D compressible Navier-Stokes Solver

<p align=center>
<img src=https://github.com/sscollis/lns3d/blob/master/docs/vortex-quarter.gif>
<br>Compressible vortex rebound from a solid wall</p>

## Background

Solves the compressible Navier-Stokes equations with 3-components (u,v,w)
and 2.5 dimensions (x,y,z) were the solution must be a single Fourier mode.

LNS3d uses a fourth-order finite-difference method and is generally designed 
to resolve highly sensitive flow phenomena such as aero-acoustics, receptivity
and linear instability.

LNS3D can also solve both linear and nonlinear 2.5d-3c problems depending on the 
problem setup but the code is specifically designed to solve for perturbations
on a baseflow to increase the fidelity of the solution.

## Building

To build you need to link (and possibly edit) the most appropriate `*.mak`
file 

    ln -s gcc.mak Makefile
    
And then

    make clean && make

## Index

Directory  |  Description
-----------|-----------------------------------------------------------
src        |  Primary source for lns3d flow solver
ssd        |  Original Cray SSD enabled source for lns3d flow solver
pre        |  Preprocessor that makes consistent metrics from a grid file
util       |  Various utilities for preparing and analyzing lns3d runs
mesh       |  Mesh generators for specific types of geometries
