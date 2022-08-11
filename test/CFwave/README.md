# Crossflow test case from Collis' PhD thesis, Ch. 4

Here we summarize the parameters for both temporal and spatial crossflow instabilty in a 
mildly compressible $M=0.3$ Falkner-Skan-Cooke boundary layer swept at $\theta=45^\circ$.  
Additional parameters are:
```bash
M = 0.3, Re = 400, Pr = 1.0
theta = 45
beta_h = 1
T_w = T_0 
```
Temporal and Spatial parameters are summarized below.  Please see the individual `README.md` in each directory 
for more details.

## Temporal case:
```bash
k_x = -0.287436451
k_z = 0.35
\omega_{lst} = 0.006533585 i
\Delta t = 0.047643

x_{min} = 0, x_{max} = 21.85938938962055
y_{min} = 0, y_{max} = 80

N_x = 20 N_y = 127
N_x = 80 N_y = 127
```

## Spatial case:
```bash
\omega = 0
k_z = 0.35
k_x_{lst} = -2.8831962907615E-001 -1.3854663677328E-002 i

x_{min} = 0, x_{max} = 217.92429906047406
y_{min} = 0, y_{max} = 80

N_x = 201 N_y = 127
```

S. Scott Collis\
flow.physics.simulation@gmail.com
