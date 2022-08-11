# Acoustic Scattering from Circular Cylinder

This is from Collis PhD Thesis, Ch. 4.

<p align=center>
<img src=https://github.com/sscollis/lns3d/blob/master/test/cyl/v.png>
<br>Acoustic wave scattered off a circular cylinder, contours of vertical 
velocity.</p>

## Run

To run, simply use the script (checking the paths to various executables).

```bash
./run.sh
```

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
4. Now trying with just one iteration.  This holds on for 4000 steps
   but it is possible that the instability could arise.

S. Scott Collis\
flow.physics.simulation@gmail.com
