# This program is developped for academic purposes 
# GNU General Public License
# copyright: jorge.brandle@coria.fr

import numpy as np
from Output import *


# Physical parameters
L = 200.0 # m
u = 1 # m/s
tmax = 1000 # s

# Boundary conditions
c_g = 0.0

# Numerical parameters
Np = 300
CFL = 0.8

# Output parameters
dtAff = 200.0 # s

# Initialise other parameters
dx = L/Np # m
dt = CFL * dx / u # s
c = np.zeros(Np)
F = np.zeros(Np)
x = np.linspace(0,L,Np)
c[0] = c_g
t = 0

# Initial condition
# Fonction Porte
c[1:Np] = (x[1:Np]<L/2+L/10)*(x[1:Np]>L/2-L/10)

# For animated figures
savefigC(c,t,x)

# Main code
while t<tmax:
    # Update fields
    # RK2 Step 1
    F = u * c
    cdemi = np.copy(c) 
    cdemi[1:Np] = cdemi[1:Np] + 0.5*dt/dx*(F[0:Np-1] - F[1:Np])
    # Boundary condition 
    # cdemi[0] = c_g   # Diriclet
    cdemi[0] = cdemi[Np-1] # Periodique
    # RK2 Step2
    F = u * cdemi
    c[1:Np] = c[1:Np] + dt/dx*(F[0:Np-1] - F[1:Np])
    # Boundary condition 
    # c[0] = c_g   # Diriclet
    c[0] = c[Np-1] # Periodique
    
    t = t + dt
    
    # Plot
    if (abs(t % dtAff) < dt):
        savefigC(c,t,x)
