# This program is developped for academic purposes 
# GNU General Public License
# copyright: jorge.brandle@coria.fr

from matplotlib.tri import Triangulation
from MeshConnectivities import *
import numpy as np

class Mesh:
     def __init__(self, xnodes, ynodes):
        self.Nrefin = 3
        self.triang = Triangulation(xnodes, ynodes)
        self.Area = ComputeArea(self.triang)
        self.xcellcen, self.ycellcen, self.triangcellcent = Compute_CellCen(self.triang)
        self.nx, self.ny, self.dS = ComputeNormals(self.triang)
        self.Ltri, self.Rtri, self.indEdge, self.direc = FindDirections(self.triang)
        self.xcellcen_fine, self.ycellcen_fine, self.triang_fine, self.father, self.triangcellcent_fine = ComputeFine(self.triang,self.Nrefin)
        self.xedgecen, self.yedgecen = Compute_EdgeCen(self.triang)




def InitSquare(Nx,Ny,Lx,Ly):
    x = np.linspace(0,Lx,Nx) - Lx/2.0
    y = np.linspace(0,Ly,Ny) - Ly/2.0
    X, Y = np.meshgrid(x,y)
    X = X #+ 0.1*np.random.random(np.shape(X))*Lx/Nx
    Y = Y #+ 0.1*np.random.random(np.shape(Y))*Ly/Ny
    xnodes = X.reshape((Nx*Ny))
    ynodes = Y.reshape((Nx*Ny))
    return Mesh(xnodes, ynodes)


def InitCircle(n_angles,n_radii,rmin,rmax):
    radii = np.linspace(0, rmax, n_radii)
    angles = np.linspace(0, 2*np.pi, n_angles, endpoint=False)
    angles = np.repeat(angles[..., np.newaxis], n_radii, axis=1)
    angles[:, 1::2] += np.pi/n_angles
    x = (radii*np.cos(angles)).flatten()
    y = (radii*np.sin(angles)).flatten()
    return Mesh(x, y)
