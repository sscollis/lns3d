# LNS3D compressible Navier-Stokes Solver

<p align=center>
<img src=https://github.com/sscollis/lns3d/blob/master/docs/vortex.gif>
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
file.  For example, using GCC Gfortran on Mac OS X (Darwin), one would use

    cd src
    ln -s gcc.mak Makefile
    
And then

    make clean && make
    
The executable is `lns3d`.

## Auxiliary Tools

The following tables describes the supporting tools for `lns3d` that are
required for complete end-to-end workflows and analysis.

Directory  |  Description
-----------|---------------------------------------------------------------
`src`      |  Primary source for lns3d flow solver
`ssd`      |  Original Cray SSD enabled source for lns3d flow solver
`nossd`    |  Original Cray version of lns3d flow solver with SSD disabled
`pre`      |  Preprocessor that makes consistent metrics from a grid file
`util`     |  Various utilities for preparing and analyzing lns3d runs
`mesh`     |  Mesh generators for specific types of geometries
`test`     |  Simple tests and examples of `lns3d` and tools

Each of these tool directories contains a separate `README` file that
describes how to build and use them.  

Finally, the `test` directory contines several examples demonstrating
workflows and capabilities of `lns3d` that can both help assess that 
your build are working correctly and show new users how to setup and
run problems.

---

S. Scott Collis\
Wed Feb 26 06:32:22 MST 2020
