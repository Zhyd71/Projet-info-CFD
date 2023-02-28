# This program is developped for academic purposes 
# GNU General Public License
# copyright: jorge.brandle@coria.fr

from matplotlib.tri import Triangulation, UniformTriRefiner,\
    CubicTriInterpolator, TriRefiner
import numpy as np
import math
from Output import *

# Physical parameters
tmax = 700 # s

# Numerical parameters
CFL = 0.4

# Output parameters
dtAff = 10.0 # s

# Make mesh
print("Making mesh")
from MeshGeneration import *
#mesh = InitSquare(15,15,100.0,100.0)
mesh = InitCircle(25,25,0,50)


# Initialise density
print("Init density")
from InitScalar import *
c = init_c(zalesak,mesh)

# Compute velocities (centered on edges)
print("Init velocity")
from InitVector import *
u, v = InitVortex(mesh, np.pi/314,0,0)

# Compute dF = n.v*dS
print("Compute dF")
def ComputeDf(u,v,mesh):
    return (u*mesh.nx + v*mesh.ny)*mesh.dS
    #dF = np.zeros(np.size(mesh.triang.edges[:,0]))
    #for i,edge in enumerate(mesh.triang.edges):
    #    dF[i] = (u[i]*mesh.nx[i] + v[i]*mesh.ny[i])*mesh.dS[i]
    #return dF
dF = ComputeDf(u,v,mesh)

def ComputeFluxUpwind(c,dF,mesh):
    return ((dF < 0)*(mesh.Ltri >= 0)*c[mesh.Ltri]/mesh.Area[mesh.Ltri] + (dF>0)*(mesh.Rtri >=0)*c[mesh.Rtri] / mesh.Area[mesh.Rtri])* dF 
    #F = np.zeros(np.size(mesh.triang.edges[:,0]))
    #for i,edge in enumerate(mesh.triang.edges):
    #    if (dF[i] < 0 and Ltri[i] >= 0):
    #        F[i] = c[Ltri[i]] * dF[i] / Area[Ltri[i]]
    #    elif (dF[i]>0 and Rtri[i] >=0):
    #        F[i] = c[Rtri[i]] * dF[i] / Area[Rtri[i]]
    #return F
            
          
# init time
t = 0
savefigC(c,t,mesh)

# Go
print("Time loop")
while t<tmax:
    print("t=",t)
    # Compute flux
    F = ComputeFluxUpwind(c,dF,mesh)
    
    # Update
    dt = 0.5
    for i,tri in enumerate(mesh.triang.triangles):
       for j,neig in enumerate(mesh.triang.neighbors[i]):
           if neig>=0:
               c[i] = c[i] + dt*F[mesh.indEdge[i,j]]*mesh.direc[i,j]
    t = t + dt
    if (abs(t % dtAff) < dt):
        savefigC(c,t,mesh)

