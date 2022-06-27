#!env python
import matplotlib.pyplot as plt
import numpy as np
wall = np.loadtxt("wall.dat.0")
plt.title('Pressure gradient')
plt.ylabel(r"$\partial p/\partial s$")
plt.xlabel(r"$x$")
plt.plot(wall[:,0],wall[:,4],'b-',label=r"$\partial p/\partial s$")
#plt.xlabel(r"$s$")
#plt.plot(wall[:,1],wall[:,4],'b-',label=r"$\partial p/\partial s$")
plt.tick_params(axis='x',direction='in')
plt.tick_params(axis='y',direction='in')
#plt.xticks(np.arange(0,18,2))
#plt.yticks(np.arange(-0.2,1.2,.2))
plt.xscale("log")
plt.axis([0.001, 5.0e2, -0.3, 0.1])
plt.legend()
plt.savefig("wall.png")
plt.show(block=False)
