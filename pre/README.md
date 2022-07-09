# Pre-processor for LNS3d

The pre-processor for LNS3d that takes an arbitrary `grid.dat` in appropriate
unformated (or binary) format (depending on the system) and generates the
`metric.dat` file using either fourth-order finite-difference or (for the
original version) B-splines.  Both versions right now only have the finite-
difference active as the B-splines relied on IMSL.  Now that `bslib` has
interfaces that match IMSL, it should be easy to get B-splines back working.

The directories herein are:

Directory  |   Description
-----------|----------------------------------------------------------------
`src_org`  |  The original version working in `JIK` ordering with B-splines
`src_ij`   |  Updated to use `KIJ` ordering but with no B-spline support
-----------|----------------------------------------------------------------

There is a symbolic link `src -> src_ij` as the version to use.

S. Scott Collis\
flow.phyaics.simulation@gmail.com
