# This program is developped for academic purposes 
# GNU General Public License
# copyright: jorge.brandle@coria.fr

import numpy as np
from Output import *


# Physical parameters
L = 200.0 # m
tmax = 700 # s

# Boundary conditions
Dirichlet = False
Cyclic = True
if Dirichlet:
    u_g = 0.0
    u_d = 0.0

# Numerical parameters
Np = 300
CFL = 0.8

# Output parameters
dtAff = 5.0 # s

# Initialise other parameters
dx = L/Np # m
u = np.zeros(Np)
F = np.zeros(Np+1)
xn = np.linspace(dx/2,L-dx/2,Np)
t = 0

# Def fluxes 
def flux_upwind_conservative_naive(U,V):
    uedge = (U + V)/2.0
    return (uedge > 0) * uedge * U + (uedge < 0) * uedge * V 

def flux_upwind_conservative(U,V):
    uedge = (U + V)/2.0
    return (uedge > 0) * U**2 + (uedge < 0) * V**2
   
def flux_Godunov(U,V):
    ustar = (U>=V)*(((U+V)>0)*U + ((U+V)<=0)*V) + (U<V)*((U>0)*U+(V<0)*V)
    return 0.5*ustar**2

def flux(U,V):
    #return flux_Godunov(U,V)
    return flux_upwind_conservative(U,V)
    #return flux_upwind_conservative_naive(U,V)

# Initial condition
# Fonction gaussienne
u = np.exp(-((xn-L/2)/(0.1*L))**2)+0.5

# For animated figures
savefigC(u,t,xn,(np.min(u),np.max(u)))

# Main code
while t<tmax:
    # Update fields
    # Godunov
    F[1:Np] = flux(u[0:Np-1],u[1:Np])
    
    if Dirichlet: 
        F[0] = flux(u_g,u[0])
        F[Np] = flux(u[Np-1],u_d)
    elif Cyclic:
        F[Np] = flux(u[Np-1],u[0])
        F[0] = F[Np]
    
    dt = CFL * dx / np.max(u) # s
    u = u + dt/dx*(F[0:Np] - F[1:Np+1])

    t = t + dt
    
    # Plot
    if (abs(t % dtAff) < dt):
        savefigC(u,t,xn,(np.min(u),np.max(u)))
