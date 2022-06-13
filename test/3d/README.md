## Inviscid oblique acoustic wave on a square, periodic box.

<p align=center>
<img src=https://github.com/sscollis/lns3d/blob/master/test/3d/3d.png>
<br>Oblique acoustic wave, contours of density.</p>

# Mesh

Nx = 64
Ny = 63

Lx = 1
Ly = 1
Lz = 1

# Flow 

Ma = 1, Re = 0, Pr = 1

Zero mean flow.

Input      | Description
-----------|-------------------------------
`ic.inp`   | mkdist3d input 3d oblique wave
`dist.inp` | mkdist3d input Old wave
`rk4.inp`  | LNS3d input using Rk4
`cfl.inp`  | LNS3d input implicit
`grid.inp` | genmesh input file
`mean.inp` | mkmesh input file

S. Scott Collis\
flow.physics.simulation@gmail.com
