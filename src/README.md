# LNS3D 

3C-2.5D unsteady compressible, linearized Navier--Stokes code 

## Change log

### 8-28-96

To compile on Unix systems type:  `make`

### 6-25-98

I have brought this up-to-date (almost) with the SSD version on the Cray C90.
Routines that have not been fully updated will generate an error upon entry.

These include:

    impdrv3d.f90

### 7-24-98  

Changed module names stuff -> global

### 12-18-00

I have switched almost all of the routines to `(ndof,i,j)` ordering which
works much better on cache-based architecture.

Routines that still need attention:

    traces.f90
    wallbc.f90
    imprk.f90
    Lele & Poinsoit BC's

The "perl" directory contains Perl scripts that helped to quickly change the
data-structures throughout the code.

### 12-19-00

1) Updated module structure
2) Fixed bugs in rkbc and Reimann
3) Need to change the way metrics are stored
4) Need to update `traces.f90` and `wallbc.f90`
5) Need to update `grad.f90` and `grad2.f90` with versions that support full
   symmetry conditions.
