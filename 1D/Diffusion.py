import numpy as np
from Output import *

# Physical parameters
L = 100.0 # m
k = 1 # 
tmax = 70 # s

# Boundary conditions
T_g = 0.0
T_d = 0.0

# Numerical parameters
Np = 300
CFL = 0.45

# Initialise other parameters
dx = L/Np # m
dt = CFL * dx*dx / k # s
T = np.zeros(Np)
F = np.zeros(Np-1)
x = np.linspace(0,L,Np)
T[0] = T_g
T[Np-1] = T_d
t = 0

# Initial condition
# T[1:Np] = (x[1:Np]<(L/2+L/30))*(x[1:Np]>L/2-L/30) 
for i in range(Np):
    if (x[i] < L/2+L/30) and (x[i]>L/2-L/30):
        T[i] = 1.
    else:
        T[i] = 0. 

savefigC(T,t,x)

# Main code
while t<tmax:
    # Update fields
    F = k * (T[1:Np]-T[0:Np-1])/dx
    T[1:Np-1] = T[1:Np-1] + dt/dx*(F[1:Np-1] - F[0:Np-2])
    t = t + dt
    savefigC(T,t,x)
    
