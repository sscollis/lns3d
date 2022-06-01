# Spatial Crossflow Vortex 

<p align=center>
<img src=https://github.com/sscollis/lns3d/blob/master/test/CFwave/spatial/v.png>
<br>Compressible temporally growing stationary crossflow vortex, contours 
of vertical velocity.</p>

## Problem Setup

From Collis, PhD Thesis, $\S4.2.2$, sweep angle $\theta=45^\circ$, 
Hartree parameter $\beta_h = 1$, $M = 0.3$, $Re = 400$, $Pr = 1$ 
with $T_w = T_0$. We use a stationary $\omega=0$ crossflow mode with
and $k_z = 0.35$ with LST yeilding $k_{x_{lst}} = -0.288319629 - 0.013854663 i$
whic means the spatial growth rate is $\sigma_{lst} = 0.013854663$.

```bash
Re =  4.000000E+02, Ma =  3.000000E-01, Pr =  1.000000E+00
Omega = ( 0.0000000000000E+000, 0.0000000000000E+000)
Alpha = (-2.8831962907615E-001,-1.3854663677328E-002)
Beta  = ( 3.5000000000000E-001, 0.0000000000000E+000)
```

The wavelength in the streamwise direction is $\lambda_x = 21.792429906047406$ 
and in the crossflow direction $\lambda_z = 17.95195802051310$

S. Scott Collis\
sscollis@gmail.com 
