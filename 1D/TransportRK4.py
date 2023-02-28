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
cTemp = np.copy(c)

# For animated figures
savefigC(c,t,x)

# Main code
while t<tmax:
    # Update fields
    # RK4 - Step1
    F1 = u * c
    cTemp[1:Np] = c[1:Np] + 0.5*dt/dx*(F1[0:Np-1] - F1[1:Np])
    cTemp[0] = cTemp[Np-1]
    F2 = u * cTemp
    cTemp[1:Np] = c[1:Np] + 0.5*dt/dx*(F2[0:Np-1] - F2[1:Np])
    cTemp[0] = cTemp[Np-1]
    F3 = u * cTemp
    cTemp[1:Np] = c[1:Np] + dt/dx*(F3[0:Np-1] - F3[1:Np])
    cTemp[0] = cTemp[Np-1]
    F4 = u * cTemp
    F = 1./6*(F1+2*F2+2*F3+F4)
    c[1:Np] = c[1:Np] + dt/dx*(F[0:Np-1] - F[1:Np])
    c[0] = c[Np-1]

    
    t = t + dt
    
    # Plot
    if (abs(t % dtAff) < dt):
        savefigC(c,t,x)
