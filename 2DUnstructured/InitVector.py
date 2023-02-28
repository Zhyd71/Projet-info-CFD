# This program is developped for academic purposes 
# GNU General Public License
# copyright: jorge.brandle@coria.fr

import numpy as np
def InitVortex(mesh,omega,x0,y0):
    u = np.zeros(np.size(mesh.triang.edges[:,0]))
    v = np.zeros(np.size(mesh.triang.edges[:,0]))
    for i,edge in enumerate(mesh.triang.edges):
        u[i] = omega*(mesh.yedgecen[i]-y0)
        v[i] = omega*(x0-mesh.xedgecen[i])
    return u, v
