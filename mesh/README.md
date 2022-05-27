## LNS Mesh tools

The following mesh generators all function however, they use the old JI 
metric format so you have to use the -ms flag on updated codes.

Code       | Description
-----------|--------------------------------------------------------
`level`    | make level-set grid given body.dat
`cyl`      | make level-set mesh for circular cylinder (IJ format)
`ogrid`    | make an ogrid for a cylinder or airfoil
`confpc`   | make a conformal grid for parabolic cylinder
`pc`       | make a body-fitted mesh for parabolic cylinder
`mse`      | modified super-ellipse mesh
`circ`     | write out x, y for a circle
`cinterpc` | interpolates from one mesh to another with same mapping
`interpc`  | interpolates from one mesh to another

## Building

The basic codes are built simply by

    ln -s gcc.mak Makefile
    make all

and this works on platforms (including MacOS) with gfortran installed.

## Building with Numerical Recipes 

Note that there are dependencies on commercial code that are, of
course, not distributed here so to build you either need a 
license to those codes or will need to replace
them with equivalents.

To use the commercial (Numerical Recipes) code, build using:

    ln -s gcc.mak Makefile
    env USE_NR=1 LIBNR_DIR=$HOME/git/NR-utilities make all

where `LIBNR_DIR` points to the directory where you have `libnr.a` containing
Numerical-Recipes routines.

This requires that you have two additional source files (not 
distributed here) `nr_odeint.f` and `nr_spline.f`. Both would be
easy to replace with public domain software and this is 
encouraged as a future update.

S. Scott Collis\
sscollis@gmail.com
