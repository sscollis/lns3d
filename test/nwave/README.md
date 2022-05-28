# Simple "nonlinear" acoustic wave test

This is a "nonlinear" acoustic wave test.  The domain is from $x \in [0,500]$; 
$y \in [0:100]$ and the initial condition is a Guassian right running 
acoustic pulse positioned at $x = 250$ with $\sigma = 25$ and $amp=1.0e-5$.  
It makes a useful test of the nonlinear solver.

I have verified that this works with all solvers:
```bash
lns3d
cns2d_prim
cns2d_cons
```

S. Scott Collis\
scollis@gmail.com
