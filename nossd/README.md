# Original Reference Version of LNS3D

This version does not use the Cray SSD but otherwise can run on 
current CPUs.

Do use with caution as it is not routinely tested.

---
8-28-96

LNS3D: 3D unsteady compressible, linearized Navier--Stokes code 

To compile on Unix systems type:  make

6-25-98

I have brought this up-to-date (almost) with the SSD version on the Cray C90.
Routines that have not been fully updated will generate an error upon entry.

These include:

    impdrv3d.f90

7-24-98  

Changed module names stuff -> global

3-14-2020

This builds and runs on Macs with Gfortran.  While it is slow (due to 
cache misses and no OpenMP parallelism, it is the only version that
is consistent and fully up-to-date relative to physics and 
boundary-conditions.
