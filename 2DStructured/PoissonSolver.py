# This program is developped for academic purposes 
# GNU General Public License
# copyright: jorge.brandle@coria.fr
import numpy as np
import scipy.sparse as sp 
import scipy.sparse.linalg  
import scipy.sparse.linalg as spla

def CreatePoissonMatrix(mesh,BC):
     # Right hand side
     b = np.zeros(mesh.Nx*mesh.Ny)
     # Diagonals :
     Diag =-2.*(1./mesh.dx**2+1./mesh.dy**2) * np.ones(mesh.Nx*mesh.Ny)
     Diagxp = 1./mesh.dx**2 * np.ones(mesh.Nx*mesh.Ny-1)
     Diagxm = 1./mesh.dx**2 * np.ones(mesh.Nx*mesh.Ny-1)
     Diagyp = 1./mesh.dy**2 * np.ones(mesh.Nx*mesh.Ny-mesh.Nx)
     Diagym = 1./mesh.dy**2 * np.ones(mesh.Nx*mesh.Ny-mesh.Nx)
     for i in range(mesh.Ny-1):
        Diagxp[(i+1)*mesh.Nx-1] = 0
        Diagxm[(i+1)*mesh.Nx-1] = 0
     # Boundary conditions for p
     if BC["type"]["Left"] in ["Wall", "MovingWall", "Inflow"]: # Neumann
         for i in range(mesh.Ny):
             Diag[i*mesh.Nx] += 1./mesh.dx**2
     elif BC["type"]["Left"] == "Outflow":
         for i in range(mesh.Ny):
             b[i*mesh.Nx] = -1./mesh.dx**2 

     if BC["type"]["Right"] in ["Wall", "MovingWall", "Inflow"]: # Neumann
         for i in range(mesh.Ny):
             Diag[(i+1)*mesh.Nx-1] += 1./mesh.dx**2
     elif BC["type"]["Right"] == "Outflow":
         for i in range(mesh.Ny):
             b[(i+1)*mesh.Nx-1] = -1.0/mesh.dx**2

     if BC["type"]["Bottom"] in ["Wall", "MovingWall", "Inflow"]: # Neumann
         for i in range(mesh.Nx):
             Diag[i] += 1./mesh.dy**2
     elif BC["type"]["Bottom"] == "Outflow":
         for i in range(mesh.Nx):
             b[i] = -1./mesh.dy**2

     if BC["type"]["Top"] in ["Wall", "MovingWall", "Inflow"]: # Neumann
         for i in range(mesh.Nx):
             Diag[mesh.Ny*mesh.Nx-i-1] +=  1./mesh.dy**2
     elif BC["type"]["Top"] == "Outflow":
         for i in range(mesh.Nx):
             b[mesh.Ny*mesh.Nx-i-1] = -1./mesh.dy**2

     Diagonals = [Diag, Diagxp, Diagxm, Diagyp, Diagym]
     A = sp.diags(np.asarray(Diagonals), [0,1,-1,mesh.Nx, -mesh.Nx]).tocsc()
     if (mesh.Nx*mesh.Nx<500):
         import matplotlib.pyplot as plt
         plt.matshow(A.todense())
         plt.show()
         plt.matshow(np.reshape(b,(mesh.Nx,mesh.Ny),order="F"))
         plt.show()
         plt.savefig("Solution/PoissonMatrix.png")
         plt.close()
     return A, b   

def PoissonSolver(A,b,bp,mesh,nit=100,LU=None):
    b = np.reshape(b,(mesh.Nx*mesh.Ny),order="F")+bp
    # x = sp.linalg.spsolve(A,b)
    # x, info = sp.linalg.bicg(A,b)
    # x, info = sp.linalg.lgmres(A,b)
    # print np.max(np.abs(A*x-b))
    x = LU(b)    
    return np.reshape(x,(mesh.Nx,mesh.Ny),order="F")
    
def Factorize(A):
    return spla.factorized(A) # Solve only once

#def PoissonSolver(b,BC,mesh,nit=100):
#    p = np.zeros_like(mesh.Xn)
#    pn = np.zeros_like(mesh.Xn)
#
#
#    # Jacobi iterative method
#    for q in range(nit):
#        pn = p.copy()
#        p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * mesh.dx**2 + 
#                          (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * mesh.dy**2) /
#                          (2 * (mesh.dx**2 + mesh.dy**2)) -
#                          mesh.dx**2 * mesh.dy**2 / (2 * (mesh.dx**2 + mesh.dy**2)) * 
#                          b[1:-1,1:-1])      
#
#        
#        p[:, -1] = p[:, -2] ##dp/dy = 0 at BC top
#        p[0, :] = p[1, :]   ##dp/dy = 0 at BC left
#        p[:, 0] = p[:, 1]   ##dp/dx = 0 at BC bottom
#        p[-1, :] = 0        ## p = 0 at y = 2
#    
#    return p-p[0,0]
        
