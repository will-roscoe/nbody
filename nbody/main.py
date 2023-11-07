import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mpl
fig = plt.figure()
ax3d = fig.add_subplot(111, projection='3d')
pos, rad = (2,4,-5), 4.56
ax3d.plot_wireframe()

import numpy as np
def sphere(pos, radius):
  for c, r in zip(pos, radius):
    # draw sphere
    u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:50j]
    x = r*np.cos(u)*np.sin(v) -c[0]
    y = r*np.sin(u)*np.sin(v) - c[1]
    z = r*np.cos(v) -c[2]
    return x,y,z

 ### ^^ returns the arrays.
