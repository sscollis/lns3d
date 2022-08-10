# Original SSD Reference Version of LNS3D

This version will error out due to Cray SSD but otherwise can build on
current CPUs.

Do use with caution as it is not routinely tested -- see the 'nossd'
directory for a reference version that runs.

Note:  
  1.  The input specification has changed for the current version in ‘src‘ so beware incompatibilities.

S. Scott Collis\
flow.physics.simulation@gmail.com

---

8-28-96

LNS3D: 3D unsteady compressible, linearized Navier--Stokes code 

To compile on Unix systems type:  make

This is the original CRAY SSD version using JI ordering.

3-14-2020

Note that it will build on Unix systems with Gfortran but it
will error out with a no SSD present message
