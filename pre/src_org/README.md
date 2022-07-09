# LNS3D pre-processor (old JI version)

This is the pre-processor for LNS calculations that uses the deprecated 
`JI` ordering back for Cray vector processors.  Takes an arbitrary grid
file as input and generates the metrics using a finite difference scheme
or B-splines.

Now use the version in the `src_ij` directory but it only supports finite-
difference right now.   

TODO:  add back in the B-spline approach to `src_ij` which should be easy 
as I have duplicated the IMSL B-spline routines using `slatec`...

S. Scott Collis\
flow.physics.simulation@gmail.com
