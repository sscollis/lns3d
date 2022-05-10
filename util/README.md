# LNS utilities

This directory contains the utility programs that support LNS.

Note that with few exceptions, it is assumed that these codes are compiled 
using the -r8 (or equivalent) option in the compiler to cast all reals to
doubles.  This is supported in the SGI makefile (sgi.mak)

## Building



Notation | Description
---------|---------------------------------------------------------------
IO       | input/output has switch to IJ but internal is still JI
NO       | no changes have been made
YES      | fully IJ throughout
N/A      | Not applicable

Executable  |   ndof,i,j   |  Comment 
------------|--------------|--------------------------------------------
conv.sgi    |    N/A       |  Convert Plot3D format from SGI double fortran unformated to single binary
csubwave    |    NO        |  Subtract off a wave from a complex field
dirp3d      |    N/A       |  Convert HDIR field to Plot3d format
genmesh     |    YES       |  Make simple rectangular mesh
initial     |    NO        |  Make simple initial conditions for LNS
inter       |    NO        |  B-spline interpolation of a LNS field from one mesh to another in computational space.
lpost       |    NO        |  LNS post processor for linear calculations
lpost3d     |    YES       |  LNS post processor for complex calculations
mkamp       |    NO        |  Determine amplitude of forced acoustic wave
mkdist      |    IO        |  LNS disturbance initial condition for 2d flows
mkdist3d    |    YES       |  LNS disturvance initial condition for 3d flows
mkini       |    YES       |  LNS initial condition for 2d flows
mkmean      |    IO        |  Make simple mean flows for LNS
mkpot       |    NO        |  Extract potential boundary values from an LNS restart file
mkvortex    |    YES       |  Make a compressible Oseen vortex IC
nconvert    |    NO        |  Put an eigenfunction on a grid
npost       |    YES       |  LNS post processor for nonlinear calculations. This utility can also switch indices on metrics and convert from/to conservative variables
p3dlns3d    |    NO        |  Convert a Plot3d file into LNS field
r4tor8      |    N/A       |  Convert a r*8 Plot3d grid file to r*4
spost       |    NO        |  LNS post processor for swept, nonlinear, calcs
subwave     |    NO        |  Remove a forcing wave from LNS field
ij2ji       |    N/A       |  Read LNS field in IJ and write in JI
ji2ij       |    N/A       |  Read LNS field in JI and write in IJ

Scripts     |   ndof,i,j   |    Comment 
------------|--------------|--------------------------------------------
conv        |   N/A        |    Runs conv.sgi for multiple fields
gconv       |   N/A        |    Runs conv.sgi -g for multiple grids
ppost       |   N/A        |    Runs npost and preplot to convert and LNS field file to a tecplot.plt file
movie_lns3d |   N/A        |    Tecplot macro for LNS3D movies
movie_cns2d |   N/A        |    Tecplot macro for CNS2D movies


Subroutines |   ndof,i,j   |    Comment 
------------|--------------|--------------------------------------------
cfilter     |   N/A        |    Filter a complex field in one direction
error       |   N/A        |    General error handler
filter      |   N/A        |    Filter a real field in one direction
grad        |   YES        |    Take spatial gradient of a field
grad2       |   YES        |    Take second derivative of field in space
wdata       |   N/A        |    Write binary Plot3d solution file
wgrid       |   N/A        |    Write binary Plot3d grid file

Similar to the lns3d mesh generators, there are minor dependencies on
commercial math libraries that the user will need to supply to successfully build
