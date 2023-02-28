# This program is developped for academic purposes 
# GNU General Public License
# copyright: jorge.brandle@coria.fr

import numpy as np
from Output import *

# Physical parameters
tmax = 628 # s

# Boundary conditions
Dirichlet = False
Neumann = False
Cyclic = True

if Dirichlet: 
    c_g = 0.0
    c_d = 0.0
    c_h = 0.0
    c_b = 0.0

# Numerical parameters

CFL = 0.4

# Make mesh
from MeshGeneration import *
mesh = Mesh(100.0,100.0,300,300)

# Output parameters
dtAff = 30. # s

# Initialise other parameters

Fx = np.zeros_like(mesh.Xex)
Fy = np.zeros_like(mesh.Xey)
u = np.zeros_like(mesh.Xex)
v = np.zeros_like(mesh.Xey)

t = 0

# Initial condition 
from InitScalar import *
c = init_c(zalesak, mesh) # Zalesak Disk

u =   np.pi/314*mesh.Xex
v =  -np.pi/314*mesh.Yey

savefigC(c,t,mesh)

# Main code
while t<tmax:
    # Compute flux
    Fx[1:mesh.Nx,:] = (u[1:mesh.Nx,:]>0)*u[1:mesh.Nx,:]*c[0:mesh.Nx-1,:] + (u[1:mesh.Nx,:]<0)*u[1:mesh.Nx,:]*c[1:mesh.Nx,:]
    Fy[:,1:mesh.Ny] = (v[:,1:mesh.Ny]>0)*v[:,1:mesh.Ny]*c[:,0:mesh.Ny-1] + (v[:,1:mesh.Ny]<0)*v[:,1:mesh.Ny]*c[:,1:mesh.Ny]
   
    # Boundary conditions
    if Dirichlet: 
        Fx[0,:] = (u[0,:]>0)*u[0,:]*c_g + (u[0,:]<0)*u[0,:]*c[0,:]
        Fx[mesh.Nx,:] = (u[mesh.Nx,:]>0)*u[mesh.Nx,:]*c[mesh.Nx-1,:] + (u[mesh.Nx,:]<0)*u[mesh.Nx,:]*c_d
        Fy[:,0] = (v[:,0]>0)*v[:,0]*c_h + (v[:,0]<0)*v[:,0]*c[:,0]
        Fy[:,mesh.Ny] = (v[:,mesh.Ny]>0)*v[:,mesh.Ny]*c[:,mesh.Ny-1] + (v[:,mesh.Ny]<0)*v[:,mesh.Ny]*c_b
    elif Neumann:
        Fx[0,:] = Fx[1,:]
        Fx[mesh.Nx,:] = Fx[mesh.Nx-1,:]
        Fy[:,0] = Fy[:,1]
        Fy[:,mesh.Ny] = Fy[:,mesh.Ny-1]
    elif Cyclic:
        Fx[0,:] = Fx[mesh.Nx-1,:]
        Fx[mesh.Nx,:] = Fx[1,:]
        Fy[:,0] = Fy[:,mesh.Ny-1]
        Fy[:,mesh.Ny] = Fy[:,1]        

    
    dt = CFL * np.min(mesh.dmin/max(np.max(np.abs(u)),np.max(np.abs(v)))) # s
    # Update fields
    c = c + dt/mesh.dx*(Fx[0:mesh.Nx,:] - Fx[1:mesh.Nx+1,:])+dt/mesh.dy*(Fy[:,0:mesh.Ny]-Fy[:,1:mesh.Ny+1])
    t = t + dt
    # savefig
    if (abs(t % dtAff) < dt):
        savefigC(c,t,mesh)

savefigC(c,t,mesh)
