# LNS Grid tools

The following grid generators all function however, they use the old JI 
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

Now working on Mac with gfortran.  

Note that there are dependencies on commercial code that are, of
course, not distributed here so to build you either need a 
license to those codes or will need to replace
them with equivalents.

S. Scott Collis\
Sat Feb 15 16:32:08 MST 2020
