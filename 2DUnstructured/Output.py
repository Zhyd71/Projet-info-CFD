# This program is developped for academic purposes 
# GNU General Public License
# copyright: jorge.brandle@coria.fr

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

def savefigC(c,t,mesh):
    # Make finer mesh to plot centered volumes
    c_fine = 1.0*np.zeros_like(mesh.triang_fine.triangles[:,0])
    for i, tri in enumerate(mesh.triang_fine.triangles):
        c_fine[i] =  c[mesh.father[i]]
    # Plot
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    ax.triplot(mesh.triang, color='0.8',  alpha=.5)
    cmap = cm.get_cmap(name='bwr', lut=None)
    levels = np.linspace(np.min(c_fine), np.max(c_fine), num=100)
    tricont = plt.tricontourf(mesh.triangcellcent_fine, c_fine, cmap=cmap, levels=levels, extend='both')
    cbar = fig.colorbar(tricont)
    plt.savefig("Solution/c"+'{:07.3f}'.format(t)+".png")
    plt.close()

def savefigV(u,v,t,mesh):
    # Plot
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    ax.triplot(mesh.triang, color='0.8',  alpha=.5)
    Max = np.max(np.sqrt(u**2+v**2))
    plt.quiver(mesh.xedgecen, mesh.yedgecen,u/Max,v/Max, scale=50,scale_units='width')
    plt.savefig("Solution/V"+'{:07.3f}'.format(t)+".png")
    plt.close()
