# LNS3D compressible linearized Navier-Stokes Solver

## Background

Solves the compressible linearized Navier-Stokes equations with 3-components (u,v,w)
and 2 dimensions (x,y).  This is sometimes refered to as 2d-3c.

LNS3d uses a fourth-order finite-difference method and is generally designed 
to resolve highly sensitive flow phenomena such as aero-acoustics, receptivity
and linear instability.

## Index

Directory  |  Description
-----------|-----------------------------------------------------------
src        |  Primary source for lns3d flow solver
ssd        |  Original Cray SSD enabled source for lns3d flow solver
pre        |  Preprocessor that makes consistent metrics from a grid file
util       |  Various utilities for preparing and analyzing lns3d runs
mesh       |  Mesh generators for specific types of geometries
