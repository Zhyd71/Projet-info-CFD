# This programm is developped for academic purposes 
# GNU General Public License
# copyright: jorge.brandle@coria.fr

import numpy as np
from Output import *
from PoissonSolver import *

# Physical parameters
Re = 10.

# Numerical parameters
CFL = 0.9
dt = 0.001/4
t = 0.0
tmax = 1.5
# Make mesh 
from MeshGeneration import *
mesh = Mesh(2.0,1.0,50,25)

# Output parameters 
dtAff = 0.01 # s

# Boundary condition 
# Cavity
#BC = {"type":{"Left": "Wall", "Right": "Wall","Bottom": "Wall", "Top": "MovingWall"},
#      "u":{"Top": 1.0}, 
#      "v":{}}
# Forcing Vortex
# BC = {"type":{"Left": "MovingWall", "Right": "MovingWall","Bottom": "MovingWall", "Top": "MovingWall"},
#      "u":{"Top": 1.0, "Bottom": -1.0},
#      "v":{"Left": 1.0, "Right": -1.0}}
# Channel
BC = {"type":{"Left": "Inflow", "Right": "Outflow","Bottom": "Wall", "Top": "Wall"},
      "u":{"Left": 1.0}, 
      "v":{"Left": 0.0},
      "p":{"Right": 0.0}}
#BC = {"type":{"Left": "Wall", "Right": "Wall","Bottom": "Outflow", "Top": "Inflow"},
#      "u":{"Top":  0.0}, 
#      "v":{"Top": -1.0}}

# Initial condition
from InitNS import *
u, v, p = InitRest(mesh)

if BC["type"]["Left"] == "Inflow":
    u[0,:] = BC["u"]["Left"]
if BC["type"]["Top"] == "Inflow":
    v[:,mesh.Ny] = BC["v"]["Top"]

# Define fluxes
from NS_fluxes import GRAD, STRESS, TRANS, DIV

# Create Matrix for Poisson Solvre
A, bp = CreatePoissonMatrix(mesh, BC)
LU = Factorize(A) # Solve only once

# Print initial condition
savefigNS(u,v,p,t,mesh)

# Start time loop
while t<tmax:
    ustar, vstar = np.copy(u), np.copy(v)
    utrans, vtrans = TRANS(u,v,mesh,BC)
    ustress, vstress = STRESS(u,v,mesh,BC)

    ustar[1:-1,:] = u[1:-1,:] + dt*(utrans+1.0/Re*ustress)
    vstar[:,1:-1] = v[:,1:-1] + dt*(vtrans+1.0/Re*vstress)
    b = DIV(ustar,vstar,mesh)
    p = PoissonSolver(A,b,bp,mesh,LU=LU)
    #p = PoissonSolver(b,BC,mesh)
    uprime, vprime = GRAD(p,mesh)
    u, v = ustar-uprime, vstar-vprime

    t = t + dt
    # Plot
    if (abs(t % dtAff) < dt):
        savefigNS(u,v,p,t,mesh)
