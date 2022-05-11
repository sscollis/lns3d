## Inviscid oblique acoustic wave on a square, periodic box.

# Mesh

Nx = 64
Ny = 63

Lx = 1
Ly = 1
Lz = 1

# Flow 

Ma = 1, Re = 0, Pr = 1

Input      | Description
-----------|-------------------------------
`ic.inp`   | mkdist3d input 3d oblique wave
`dist.inp` | mkdist3d input Old wave
`rk4.inp`  | LNS3d input using Rk4
`cfl.inp`  | LNS3d input implicit
`grid.inp` | genmesh input file
`mean.inp` | mkmesh input file
