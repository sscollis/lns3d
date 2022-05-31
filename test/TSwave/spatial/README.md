# Spatial case

<p align=center>
<img src=https://github.com/sscollis/lns3d/blob/master/test/TSwave/spatial//v.png>
<br>Compressible spatially growing Tollmein-Schlichting wave, contours 
of vertical velocity.  This particular image is taken using the default
input file `lns3d.inp` which runs a frequency-domain case with outflow 
sponge to zero disturbance.</p>

## Spatial TS results

```bash
FSC: M=0.3, lambda=0, beta_h=0, Tw/T0=1, Re_\delta_1 = 1000, Pr=1

k_x = 2.2804739411412E-001  -6.5163146827854E-003  (yi=1.0, 64)
k_x = 2.2804739411185E-001  -6.5163146952733E-003  (yi=1.0, 96)

\omega = 0.08
\lambda_{TS} = 2\pi/k_x = 27.552103068969444
```

The eigenfunction from `stab` is located in `efun.pro`

## Notes

1. The basic test is in the frequency domain with $\omega=0.08$
2. The domain is $10\lambda_{ts}$ long and the sponge is $2\lambda_{ts}$
3. As in the thesis, the initial run uses a sponge to zero disturbanced `ispg = 1`
4. You can also sponge to the eigenfunction using `ispg = 4`
5. The right BC can either be zero distrubance `right = 0` or the eigenfunction `right = 4`
6. The left BC is the prescripbed eigenfunction `left = 4`
7. You can try various combination of explicit and implicit methods
8. You can also try the time unsteady solver however the inflow disturbanced don't appear
   to be correct for that case?

S. Scott Collis\
sscollis@gmail.com
