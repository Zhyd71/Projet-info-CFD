from matplotlib.tri import Triangulation, UniformTriRefiner,\
    CubicTriInterpolator, TriRefiner
import numpy as np
import math
from Output import *

# Physical parameters
tmax = 1 # s
nu = 1.

# Numerical parameters
CFL = 0.4

# Output parameters
dtAff = 0.001 # s

# Make mesh
print("Making mesh")
from MeshGeneration import *
mesh = InitSquare(30,30,100.0,100.0)
#mesh = InitCircle(50,50,0,50)


# Initialise density
print("Init density")
from InitScalar import *
cu = init_c(cercle,mesh)
cv = init_c(cercle,mesh)

# Compute velocities (centered on edges)
print("Init velocity")
def Interp(mesh,cu,cv):
    u = np.zeros(np.size(mesh.triang.edges[:,0]))
    v = np.zeros(np.size(mesh.triang.edges[:,0]))
    for i,edge in enumerate(mesh.triang.edges):
        u[i] = (cu[mesh.Ltri[i]]+cu[mesh.Rtri[i]])/2.0
        v[i] = (cv[mesh.Ltri[i]]+cv[mesh.Rtri[i]])/2.0
    return u, v
u, v = Interp(mesh,cu,cv)

# Compute dF = n.v*dS
def ComputeDf(u,v,mesh):
    return (u*mesh.nx + v*mesh.ny)*mesh.dS
    #dF = np.zeros(np.size(mesh.triang.edges[:,0]))
    #for i,edge in enumerate(mesh.triang.edges):
    #    dF[i] = (u[i]*mesh.nx[i] + v[i]*mesh.ny[i])*mesh.dS[i]
    #return dF

def ComputeFluxUpwind(c,dF,mesh):
    return ((dF < 0)*(mesh.Ltri >= 0)*c[mesh.Ltri]/mesh.Area[mesh.Ltri] + (dF>0)*(mesh.Rtri >=0)*c[mesh.Rtri] / mesh.Area[mesh.Rtri])* dF 
    #F = np.zeros(np.size(mesh.triang.edges[:,0]))
    #for i,edge in enumerate(mesh.triang.edges):
    #    if (dF[i] < 0 and Ltri[i] >= 0):
    #        F[i] = c[Ltri[i]] * dF[i] / Area[Ltri[i]]
    #    elif (dF[i]>0 and Rtri[i] >=0):
    #        F[i] = c[Rtri[i]] * dF[i] / Area[Rtri[i]]
    #return F

def ComputeFluxDiff(c,mesh):
    return -nu*(c[mesh.Ltri]-c[mesh.Rtri])/(mesh.Area[mesh.Ltri] + mesh.Area[mesh.Rtri])      
          
# init time
t = 0
savefigC(np.sqrt(cu**2+cv**2),t,mesh)
savefigV(u,v,t,mesh)

# Go
print("Time loop")
while t<tmax:
    print("t=",t)
    # Compute flux
    dF = ComputeDf(u,v,mesh)
    Fu = ComputeFluxUpwind(cu,dF,mesh) + ComputeFluxDiff(cu,mesh)
    Fv = ComputeFluxUpwind(cv,dF,mesh) + ComputeFluxDiff(cv,mesh)
    # Update
    dt = 0.05
    for i,tri in enumerate(mesh.triang.triangles):
       for j,neig in enumerate(mesh.triang.neighbors[i]):
           if neig>=0:
               cu[i] = cu[i] + dt*Fu[mesh.indEdge[i,j]]*mesh.direc[i,j]
               cv[i] = cv[i] + dt*Fv[mesh.indEdge[i,j]]*mesh.direc[i,j]
    u, v = Interp(mesh,cu,cv)
    t = t + dt
    if (abs(t % dtAff) < dt):
        savefigC(np.sqrt(cu**2+cv**2),t,mesh)
        savefigV(u,v,t,mesh)

