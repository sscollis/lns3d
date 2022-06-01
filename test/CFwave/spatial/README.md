# Spatial Crossflow Vortex 

TODO:  update to spatial case

<p align=center>
<img src=https://github.com/sscollis/lns3d/blob/master/test/CFwave/temporal/v.png>
<br>Compressible temporally growing stationary crossflow vortex, contours 
of vertical velocity.</p>

## Problem Setup

From Collis, PhD Thesis, $\S4.2.1$, sweep angle $\theta=45^\circ$, 
Hartree parameter $\beta_h = 1$, $M = 0.3$, $Re = 400$, $Pr = 1$ 
with $T_w = T_0$. We use a stationary crossflow mode with $k_x = -0.287436451$
and $k_z = 0.35$ with LST yeilding $\omega_{lst} = 0.006533585 \iota$ which 
means the temporal growth rate is $\sigma = 0.006533585$

```bash
Re =  4.000000E+02, Ma =  3.000000E-01, Pr =  1.000000E+00
Omega = ( 6.3418480187508E-007, 6.5335847258858E-003)
Alpha = (-2.8743645100000E-001, 0.0000000000000E+000)
Beta  = ( 3.5000000000000E-001, 0.0000000000000E+000)
```

The wavelength in the streamwise direction is $\lambda_x = 21.85938938962055$ 
and in the crossflow direction $\lambda_z = 17.95195802051310$

S. Scott Collis\
sscollis@gmail.com 
