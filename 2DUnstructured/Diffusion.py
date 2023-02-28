from matplotlib.tri import Triangulation, UniformTriRefiner,\
    CubicTriInterpolator, TriRefiner
import numpy as np
import math
from Output import *

# Physical parameters
tmax = 700 # s
k = 1.0

# Numerical parameters
CFL = 0.4

# Output parameters
dtAff = 10.0 # s

# Make mesh
print("Making mesh")
from MeshGeneration import *
#mesh = InitSquare(30,30,100.0,100.0)
mesh = InitCircle(25,25,0,50)


# Initialise density
print("Init density")
from InitScalar import *
c = init_c(cercle,mesh)

def ComputeFluxDiff(c,mesh):
    return -k*(c[mesh.Ltri]-c[mesh.Rtri])/(mesh.Area[mesh.Ltri] + mesh.Area[mesh.Rtri])

            
# init time
t = 0
savefigC(c,t,mesh)

# Go
print("Time loop")
while t<tmax:
    print("t=",t)
    # Compute flux
    F = ComputeFluxDiff(c,mesh)
    
    # Update
    dt = 0.5
    for i,tri in enumerate(mesh.triang.triangles):
       for j,neig in enumerate(mesh.triang.neighbors[i]):
           if neig>=0:
               c[i] = c[i] + dt*F[mesh.indEdge[i,j]]*mesh.direc[i,j]
    t = t + dt
    if (abs(t % dtAff) < dt):
        savefigC(c,t,mesh)

