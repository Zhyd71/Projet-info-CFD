# This program is developped for academic purposes 
# GNU General Public License
# copyright: jorge.brandle@coria.fr

import numpy as np
from Output import *

# Physical parameters
tmax = 0.2 # s
gamma = 1.4 

# Numerical parameters
CFL = 0.9

# Make mesh
from MeshGeneration import *
mesh = Mesh(1.0,1.0,200,200)

# Output parameters
dtAff = 0.01 # s

# Initialise other parameters
Fx = {"rho": np.zeros_like(mesh.Xex), "rhou": np.zeros_like(mesh.Xex), "rhov": np.zeros_like(mesh.Xex), "rhoE":  np.zeros_like(mesh.Xex)}
Fy =  {"rho": np.zeros_like(mesh.Xey), "rhou": np.zeros_like(mesh.Xey), "rhov": np.zeros_like(mesh.Xey), "rhoE":  np.zeros_like(mesh.Xey)}

t = 0

# Boundary conditions
Cyclic = True
Dirichlet = False # Not yet completely implemented

# Initial condition 
from InitEuler import *
# U = CircleShock(mesh,gamma)
U = BandeShock(mesh, gamma)
# U, Ulg, Urg, Udg, Uug = Cavity(mesh) # global boundary conditions provided

# Definition of basic equations
def update_p(U):
    return (gamma-1)*(U["rhoE"]-0.5*(U["rhou"]**2+U["rhov"]**2)/U["rho"])
p = update_p(U)

def update_c(rho, p):
    return np.sqrt(gamma*p/rho)
c = update_c(U["rho"], p)
u = U["rhou"]/U["rho"]
v = U["rhov"]/U["rho"]

# Define Fluxes
def physic_flux_x(U):
    u = U["rhou"]/U["rho"]
    p = update_p(U)
    F = {}
    F["rho"] = U["rhou"]
    F["rhou"] = u*U["rhou"] + p
    F["rhov"] = u*U["rhov"]
    F["rhoE"] = u*(U["rhoE"]+p)
    return F

def Fstar(U_l,U_r,F_l,F_r,S_l,S_r):
    return (S_r*F_l-S_l*F_r+S_r*S_l*(U_r-U_l))/(S_r-S_l)

def flux_HLL_x(Ul,Ur):
    p_l = update_p(Ul)
    p_r = update_p(Ur) # Here, for simplicity we compute twice many values

    c_l = update_c(Ul["rho"], p_l)
    c_r = update_c(Ur["rho"], p_r)
    
    u_l = Ul["rhou"]/Ul["rho"]
    u_r = Ur["rhou"]/Ur["rho"]

    S_l = (u_l - c_l<0)*(u_l - c_l)
    S_r = (u_r + c_r>0)*(u_r + c_r)
    
    Fphys_l = physic_flux_x(Ul)
    Fphys_r = physic_flux_x(Ur)
    F = {}
    for key in Ur:
        F[key] = ((S_l>=0)*Fphys_l[key] + 
            (S_l<0)*(S_r>0)*Fstar(Ul[key],Ur[key],Fphys_l[key],Fphys_r[key],S_l,S_r) +
            (S_r<=0)*Fphys_r[key])

    return F

def flux_x(Ul,Ur):
    return flux_HLL_x(Ul,Ur)

def flux_y(Ud,Uu):
   Ul = {"rho": Ud["rho"], "rhou": Ud["rhov"], "rhov": Ud["rhou"],"rhoE": Ud["rhoE"]}
   Ur = {"rho": Uu["rho"], "rhou": Uu["rhov"], "rhov": Uu["rhou"],"rhoE": Uu["rhoE"]}
   Fx = flux_x(Ul,Ur)
   return {"rho": Fx["rho"], "rhou": Fx["rhov"], "rhov": Fx["rhou"],"rhoE": Fx["rhoE"]}

savefigEuler(U,t,mesh)

# Main code
while t<tmax:
    # Compute flux
    Ul = {}
    Ur = {}
    for key in U:
        Ul[key] = U[key][0:mesh.Nx-1,:]
        Ur[key] = U[key][1:mesh.Nx,:]
    Fx_unbound = flux_x(Ul,Ur)
    Ud = {}
    Uu = {}
    for key in U:
        Ud[key] = U[key][:,0:mesh.Ny-1]
        Uu[key] = U[key][:,1:mesh.Ny]
    Fy_unbound = flux_y(Ud,Uu)    
    
    # Boundary Conditions
    if Cyclic:
        for key in U:
            Ul[key] = U[key][mesh.Nx-1,:]
            Ur[key] = U[key][0,:]
            Ud[key] = U[key][:,mesh.Ny-1]
            Uu[key] = U[key][:,0]
        Flx = flux_x(Ul,Ur)
        Frx = Flx
        Fdy = flux_y(Ud,Uu)
        Fuy = Fdy
    elif Dirichlet:
        for key in U:
            Ur[key] = U[key][mesh.Nx-1,:]
            Ul[key] = U[key][0,:]
            Uu[key] = U[key][:,mesh.Ny-1]
            Ud[key] = U[key][:,0]        
        Flx = flux_x(Ulg,Ul)
        Frx = flux_x(Ur,Urg)
        Fdy = flux_y(Udg,Ud)
        Fuy = flux_y(Uu,Uug)  
  
    for key in U:
        Fx[key][1:mesh.Nx,:] = Fx_unbound[key]
        Fx[key][0,:] = Flx[key][:]
        Fx[key][mesh.Nx,:] = Frx[key][:]

        Fy[key][:,1:mesh.Ny] = Fy_unbound[key]
        Fy[key][:,0] = Fdy[key][:]
        Fy[key][:,mesh.Ny] = Fuy[key][:]    

    # CFL condition
    vmax = np.max(np.abs(u**2+v**2) + c)
    dt = CFL * mesh.dmin / vmax # s

    # Update fields
    for key in U:
        U[key] = U[key]+dt/mesh.dy*(Fy[key][:,0:mesh.Ny] - Fy[key][:,1:mesh.Ny+1])+dt/mesh.dx*(Fx[key][0:mesh.Nx,:] - Fx[key][1:mesh.Nx+1,:])

    t = t + dt
    print("t=",t)

    p = update_p(U)
    c = update_c(U["rho"], p)
    u = U["rhou"]/U["rho"]
    v = U["rhov"]/U["rho"]
    # savefig
    if (abs(t % dtAff) < dt):
        savefigEuler(U,t,mesh)

savefigEuler(U,t,mesh)
