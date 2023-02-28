# This program is developped for academic purposes 
# GNU General Public License
# copyright: jorge.brandle@coria.fr

from matplotlib.tri import Triangulation
import numpy as np

class Mesh:
     def __init__(self, Lx, Ly, Nx, Ny):
        # Grid        |-----|-----|-----| ... |-----|
        # Space     -L/2                           L/2
        # Node number    0     1     2          N-1
        # Edge number 0     1     2     3    N-1    N
        self.Nx = Nx
        self.Ny = Ny
        x = np.linspace(0,Lx,Nx+1)-Lx/2.
        y = np.linspace(0,Ly,Ny+1)-Ly/2.
        self.dx = x[1]-x[0]
        self.dy = y[1]-y[0]
        self.dmin = min(self.dx,self.dy)
        xnodes = (x[0:Nx]+x[1:Nx+1])/2.0
        ynodes = (y[0:Ny]+y[1:Ny+1])/2.0
        self.Xn, self.Yn = np.meshgrid(ynodes,xnodes)
        self.Xex, self.Yex = np.meshgrid(ynodes, x)
        self.Xey, self.Yey = np.meshgrid(y, xnodes)






