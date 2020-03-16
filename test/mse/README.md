## Modified Super Ellipse 

Mach 0.4 potential flow past a modified super ellipse using
a conformal mesh. 

### Steps to run

1. Check the paths in `run.sh`
2. Note that you need the `mse` mesh generator
3. `./run.sh`
4. `cleanup` returns the directory to original condition

Note that to enhanse convergence, this script restarts with
several different values of `sigma`.  While this works, it
is certainly not done in an optimal manner. 

### Visualization

To use `paraview` you may want to do:

    ln -s grid.dat grid.xyz

### Sample Results

#### Contours of streamwise velocity

![Streamwise velocity](https://github.com/sscollis/npot/blob/master/test/mse/u.png)

#### Computational Mesh

![Mesh](https://github.com/sscollis/npot/blob/master/test/mse/mesh.png)

#### Computational Mesh close-up

![Mesh](https://github.com/sscollis/npot/blob/master/test/mse/mesh-cu.png)

S. Scott Collis \
Wed Mar 11 06:32:22 MDT 2020
