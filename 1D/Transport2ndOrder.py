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
Dirichlet = False
Cyclic = True
if Dirichlet:
    c_g = 0.0

# Numerical parameters
Np = 300
CFL = 0.2

# Output parameters
dtAff = 0.1 # s

# Initialise other parameters
dx = L/Np # m
dt = CFL * dx / u # s
c = np.zeros(Np)
F = np.zeros(Np+1)
xn = np.linspace(dx/2,L-dx/2,Np)
t = 0

# Initial condition
# Fonction Porte
c = (xn<L/2+L/10)*(xn>L/2-L/10)

# For animated figures
savefigC(c,t,xn)

# Main code
while t<tmax:
    # Update fields
    F[2:Np+1] = u * (1.5*c[1:Np] - 0.5*c[0:Np-1])
    if Dirichlet: 
        F[1] = u * (1.5*c[0] - 0.5*c_g)
        F[0] = u * c_g
    elif Cyclic:
        F[1] = F[Np]
        F[0] = F[Np-1]
    
    c = c + dt/dx*(F[0:Np] - F[1:Np+1])

    t = t + dt
    
    # Plot
    if (abs(t % dtAff) < dt):
        savefigC(c,t,xn,(np.min(c),np.max(c)))
