# This program is developped for academic purposes 
# GNU General Public License
# copyright: jorge.brandle@coria.fr

import numpy as np
from Output import *

# Physical parameters
tmax = 628 # s

# Boundary conditions
Dirichlet = True
Cyclic = False

if Dirichlet: 
    c_g = 0.0
    c_d = 0.0
    c_h = 0.0
    c_b = 0.0

# Numerical parameters

CFL = 0.2

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

# Initial condition # Zalesak Disk
from InitScalar import *
c = init_c(zalesak, mesh)
 
u =   np.pi/314*mesh.Xex
v =  -np.pi/314*mesh.Yey

savefigC(c,t,mesh)

# Main code
while t<tmax:
    # Compute flux
    Fx[2:mesh.Nx-1,:] = (u[2:mesh.Nx-1,:]>0)*u[2:mesh.Nx-1,:]*(1.5*c[1:mesh.Nx-2,:]-0.5*c[0:mesh.Nx-3,:]) + \
                    (u[2:mesh.Nx-1,:]<0)*u[2:mesh.Nx-1,:]*(1.5*c[2:mesh.Nx-1,:]-0.5*c[3:mesh.Nx,:])
    Fy[:,2:mesh.Ny-1] = (v[:,2:mesh.Ny-1]>0)*v[:,2:mesh.Ny-1]*(1.5*c[:,1:mesh.Ny-2]-0.5*c[:,0:mesh.Ny-3]) + \
                    (v[:,2:mesh.Ny-1]<0)*v[:,2:mesh.Ny-1]*(1.5*c[:,2:mesh.Ny-1]-0.5*c[:,3:mesh.Ny])
   
    # Dirichlet Ordre 1 
    if Dirichlet:
        Fx[1,:] = (u[1,:]>0)*u[1,:]*c_g + (u[1,:]<0)*u[1,:]*c[1,:]
        Fx[mesh.Nx-1,:] = (u[mesh.Nx-1,:]>0)*u[mesh.Nx-1,:]*c[mesh.Nx-2,:] + (u[mesh.Nx-1,:]<0)*u[mesh.Nx-1,:]*c_d
        Fy[:,1] = (v[:,1]>0)*v[:,1]*c_h + (v[:,1]<0)*v[:,1]*c[:,1]
        Fy[:,mesh.Ny-1] = (v[:,mesh.Ny-1]>0)*v[:,mesh.Ny-1]*c[:,mesh.Ny-2] + (v[:,mesh.Ny-1]<0)*v[:,mesh.Ny-1]*c_b
    
        Fx[0,:] = (u[0,:]>0)*u[0,:]*c_g + (u[0,:]<0)*u[0,:]*c[0,:]
        Fx[mesh.Nx,:] = (u[mesh.Nx,:]>0)*u[mesh.Nx,:]*c[mesh.Nx-1,:] + (u[mesh.Nx,:]<0)*u[mesh.Nx,:]*c_d
        Fy[:,0] = (v[:,0]>0)*v[:,0]*c_h + (v[:,0]<0)*v[:,0]*c[:,0]
        Fy[:,mesh.Ny] = (v[:,mesh.Ny]>0)*v[:,mesh.Ny]*c[:,mesh.Ny-1] + (v[:,mesh.Ny]<0)*v[:,mesh.Ny]*c_b
    elif Cyclic:
        Fx[0,:] = Fx[mesh.Nx-3,:]
        Fx[1,:] = Fx[mesh.Nx-2,:]
        Fx[mesh.Nx-1,:] = Fx[2,:]
        Fx[mesh.Nx,:] = Fx[3,:]
        Fy[:,0] = Fy[:,mesh.Ny-3]
        Fy[:,1] = Fy[:,mesh.Ny-2]
        Fy[:,mesh.Ny-1] = Fy[:,2]        
        Fy[:,mesh.Ny] = Fy[:,3]        

    
    dt = CFL * np.min(mesh.dmin/max(np.max(np.abs(u)),np.max(np.abs(v)))) # s
    # Update fields
    c = c + dt/mesh.dx*(Fx[0:mesh.Nx,:] - Fx[1:mesh.Nx+1,:])+dt/mesh.dy*(Fy[:,0:mesh.Ny]-Fy[:,1:mesh.Ny+1])
    t = t + dt
    # savefig
    if (abs(t % dtAff) < dt):
        savefigC(c,t,mesh)

savefigC(c,t)
