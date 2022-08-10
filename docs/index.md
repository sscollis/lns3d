---
layout: default
---

# Linearized Navier-Stokes 3D (LNS3D)

LNS3D is part of the [Flow Physics & Simularion (FPS)](https://flow-physics-simulation.github.io/flow-physics-simulation/) suite of codes designed to work together to help analyize and explore fluid mechanics for both incompresible and compressible flows.  Specifically, LNS3D solves the linear (and nonlinear) compressible Navier-Stokes equations in curvalinear coordinates.   In general, linear solutions are 2d3c (harmonic in $z$) and either harmonic in $t$ or fully transient in $t$.   Nonlinear solutions are 2d in $(x,y)$ and either steady or unsteady in $t$.

LNS3D has a number of mesh generators and utilities that support the core calculations as well as the other codes in the FPS Suite.

## Types of analyses

LNS3D supports a range of linear analyses including 2d, 2d3c, and 3d 
(single mode in z) with both steady and transient capabilities. 

In addition, there is support for nonlinear analysis for 2d steady and 
unsteady flows that is used to compute base flow solutions for subsequent 
linear analysis. 

## Example: Crossflow Instabilty and Receptivity

The most thorough example of using LNS3D replicates a number of the 
calculations in Chapter 5 of Collis' PhD thesis for crossflow instabilty
and receptivity on a swept parabolic cylinder.  This geometry mimics the
leading-edge of a swept airfoil but providing an extended region of crossflow
instability growth to explore both curvature and nonparallel effects.
  1. [Compressible potential flow](https://github.com/sscollis/lns3d/blob/master/test/pcyl/sweep=35/M=0.8/README.md) at $M=0.8$
  2. [Viscous mean flow](https://github.com/sscollis/lns3d/blob/master/test/pcyl/sweep=35/M=0.8/Re=1e5/README.md) at $Re=10^5$, $Pr=1$ with sweep angle of $\theta=35^\circ$
  3. [Linear stability theory (LST)](https://github.com/sscollis/lns3d/blob/master/test/pcyl/sweep=35/M=0.8/Re=1e5/lst/README.md) at various spanwise wavenumbers with inclusion of both curvature and nonparallel effects 
  5. Comparison of direct [Linearized Navier-Stokes](https://github.com/sscollis/lns3d/blob/master/test/pcyl/sweep=35/M=0.8/Re=1e5/lns/README.md) computations with with LST where the crossflow instability eigenfunction is forced on the 
inflow to an LNS calculation.
  6. There is also a [small LNS calculation](https://github.com/sscollis/lns3d/blob/master/test/pcyl/sweep=35/M=0.8/Re=1e5/short/README.md) that is helpful
for quick experimentation, again for a forced eigenfunction on the inflow
boundary..
  7. Receptivity to a linearized Gaussian [surface roughness](https://github.com/sscollis/lns3d/blob/master/test/pcyl/sweep=35/M=0.8/Re=1e5/recep/README.md) 
is the capstone example calculation.

## Other Test Cases

The [test directory](https://github.com/sscollis/lns3d/blob/master/test) is
populated with a number of example cases (many documented in Collis' PhD
thesis) and the reader is encouraged to explore these.  Most have driver
scripts called `run.sh` and `cleanup.sh` that can be executed directly or 
with some slight modification.

## Point of contact

S. Scott Collis\
flow.physics.simulation@gmail.com
