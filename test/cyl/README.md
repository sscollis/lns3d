# Acoustic Scattering from Circular Cylinder

This is from Collis PhD Thesis, Ch. 4.

## Notes:
1. You need to take multiple passes of the implicit solver as the splitting
   error creates an instability at the leading edge
2. I have found that 3 iterations are needed for Backward Euler to be stable 
   at CFL = 2.  This seems low and has something to do with the interaction
   of the symetry BC with the wall BC at $i=1$,$j=1$.
3. Ultimately the issue is one of resolution and the use of a centered
   difference scheme.  I have upped the resolution to $n_x=512$, $n_y=256$ and
   with the backward Euler scheme with 3 iterations, this appears to be
   stable.  This is alot of iterations...

S. Scott Collis\
sscollis@gmail.com
