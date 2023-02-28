import numpy as np
from Output import *


# Physical parameters
L = 1.0 # m
tmax = 0.2 # s
gamma = 1.4 

# Boundary conditions
Cyclic = False
Dirichlet = True
 
# Numerical parameters
Np = 200
CFL = 0.8

# Output parameters
dtAff = 0.01 # s

# Initialise other parameters
dx = L/Np # m
rho = np.ones(Np)
rhou = np.zeros(Np)
rhoE = np.ones(Np)
Frho = np.zeros(Np+1)
Frhou = np.zeros(Np+1)
FrhoE =  np.zeros(Np+1)
xn = np.linspace(dx/2,L-dx/2,Np)
t = 0

# Initial condition
rho_l = 1.
rho_r = 0.125
p_l = 1.0
p_r = 0.1
u_l = 0.0
u_r = 0.0
rhou_l = rho_l*u_l
rhou_r = rho_r*u_r
rhoE_l = p_l/(gamma-1) + 0.5*rho_l*u_l**2
rhoE_r = p_r/(gamma-1) + 0.5*rho_r*u_r**2

rho[:Np/2] = rho_l
rho[Np/2:] = rho_r
rhou[:Np/2] = rhou_l
rhou[Np/2:] = rhou_r
rhoE[:Np/2] = rhoE_l
rhoE[Np/2:] = rhoE_r


# Definition of basic equations
def update_p(rho,rhou,rhoE):
    return (gamma-1)*(rhoE-0.5*rhou**2/rho)
p = update_p(rho, rhou, rhoE)

def update_c(rho, p):
    return np.sqrt(gamma*p/rho)
c = update_c(rho, p)
u = rhou/rho

# Define Fluxes
def physic_flux(rho, rhou, rhoE):
    u = rhou/rho
    p = update_p(rho,rhou,rhoE)
    Frho = rhou
    Frhou = rhou**2/rho + p
    FrhoE = (rhoE+p)*u
    return Frho, Frhou, FrhoE

def flux_HLL(rho_l,rhou_l,rhoE_l,rho_r,rhou_r,rhoE_r):
    p_l = update_p(rho_l,rhou_l,rhoE_l)
    p_r = update_p(rho_r,rhou_r,rhoE_r) # Here, for simplicity we compute twice many values

    c_l = update_c(rho_l, p_l)
    c_r = update_c(rho_r, p_r)

    def Fstar(U_l,U_r,F_l,F_r,S_l,S_r):
        return (S_r*F_l-S_l*F_r+S_r*S_l*(U_r-U_l))/(S_r-S_l)

    S_l = (rhou_l/rho_l - c_l<0)*(rhou_l/rho_l - c_l)
    S_r = (rhou_r/rho_r + c_r>0)*(rhou_r/rho_r + c_r)
    
    Frhophys_l, Frhouphys_l, FrhoEphys_l = physic_flux(rho_l, rhou_l, rhoE_l)
    Frhophys_r, Frhouphys_r, FrhoEphys_r = physic_flux(rho_r, rhou_r, rhoE_r)
    
    Frho = ((S_l>=0)*Frhophys_l + 
            (S_l<0)*(S_r>0)*Fstar(rho_l,rho_r,Frhophys_l,Frhophys_r,S_l,S_r) +
            (S_r<=0)*Frhophys_r)
    Frhou = ((S_l>=0)*Frhouphys_l + 
             (S_l<0)*(S_r>0)*Fstar(rhou_l,rhou_r,Frhouphys_l,Frhouphys_r,S_l,S_r) +
             (S_r<=0)*Frhouphys_r)
    FrhoE = ((S_l>=0)*FrhoEphys_l +
             (S_l<0)*(S_r>0)*Fstar(rhoE_l,rhoE_r,FrhoEphys_l,FrhoEphys_r,S_l,S_r) +
             (S_r<=0)*FrhoEphys_r)
    return Frho, Frhou, FrhoE 

def flux(rho_l,rhou_l,rhoE_l,rho_r,rhou_r,rhoE_r):
    return flux_HLL(rho_l,rhou_l,rhoE_l,rho_r,rhou_r,rhoE_r)

# For animated figures
savefigEuler(rho,rhou,rhoE,p,c,t,xn)

# Main code
while t<tmax:
    # Update fields
    Frho[1:Np], Frhou[1:Np], FrhoE[1:Np] = flux(rho[0:Np-1],rhou[0:Np-1],rhoE[0:Np-1],rho[1:Np],rhou[1:Np],rhoE[1:Np])

    if Cyclic:
        Frho[Np], Frhou[Np], FrhoE[Np] = flux(rho[Np-1],rhou[Np-1],rhoE[Np-1],rho[0],rhou[0],rhoE[0])
        Frho[0], Frhou[0], FrhoE[0] = Frho[Np], Frhou[Np], FrhoE[Np]
    elif Dirichlet:
       Frho[Np], Frhou[Np], FrhoE[Np] = flux(rho[Np-1],rhou[Np-1],rhoE[Np-1],rho_r,rhou_r,rhoE_r)
       Frho[0], Frhou[0], FrhoE[0] = flux(rho_l,rhou_l,rhoE_l,rho[0],rhou[0],rhoE[0])
    vmax = np.max(np.abs(u) + c)
    
    dt = CFL * dx / vmax # s
    rho = rho + dt/dx*(Frho[0:Np] - Frho[1:Np+1])
    rhou = rhou + dt/dx*(Frhou[0:Np] - Frhou[1:Np+1])
    rhoE = rhoE + dt/dx*(FrhoE[0:Np] - FrhoE[1:Np+1])

    p = update_p(rho,rhou,rhoE) 
    c = update_c(rho, p)
    u = rhou/rho

    t = t + dt
    print "t=",t

    # Plot
    if (abs(t % dtAff) < dt):
        savefigEuler(rho,rhou,rhoE,p,c,t,xn)
