# This programm is developped for academic purposes 
# GNU General Public License
# copyright: jorge.brandle@coria.fr
import numpy as np

## X-EQUATION
def DUUDX(u,mesh):
    return 0.25*((u[:-1,:]+u[1:,:])**2)

def DUVDY(u,v,mesh,BC):
    dUVdy = np.zeros((mesh.Nx-1,mesh.Ny+1))
    dUVdy[:,1:-1] = 0.25*((v[:-1,1:-1]+v[1:,1:-1])*(u[1:-1,:-1]+u[1:-1,1:]))
    # Boundary conditions
    if (BC["type"]["Bottom"] == "MovingWall"):
        dUVdy[:,0] = 0.5*((v[:-1,0]+v[1:,0])*(BC["u"]["Bottom"]))
    if (BC["type"]["Top"] == "MovingWall"):
        dUVdy[:,-1] = 0.5*((v[:-1,-1]+v[1:,-1])*(BC["u"]["Top"]))
    return dUVdy

def D2UDX2(u,mesh):
    return (u[1:,:]-u[:-1,:])/mesh.dx

def D2UDY2(u,mesh,BC):
    d2udy2 = np.zeros((mesh.Nx-1,mesh.Ny+1))
    d2udy2[:,1:-1] = (u[1:-1,1:]-u[1:-1,:-1])/mesh.dy
    # Boundary conditions
    if (BC["type"]["Bottom"] == "MovingWall"):
        d2udy2[:,0] = 2*(u[1:-1,0]-BC["u"]["Bottom"])/mesh.dy
    elif (BC["type"]["Bottom"] == "Wall"):
        d2udy2[:,0] = 2*u[1:-1,0]/mesh.dy
    if (BC["type"]["Top"] == "MovingWall"):
        d2udy2[:,-1] = 2*(BC["u"]["Top"]-u[1:-1,-1])/mesh.dy
    elif (BC["type"]["Top"] == "Wall"):
        d2udy2[:,-1] = -2*u[1:-1,-1]/mesh.dy
    return d2udy2


## Y-EQUATION
def DVVDY(v,mesh):
    return 0.25*((v[:,:-1]+v[:,1:])**2)

def DUVDX(u,v,mesh,BC):
    dUVdx = np.zeros((mesh.Nx+1,mesh.Ny-1))
    dUVdx[1:-1,:] = 0.25*((u[1:-1,:-1]+u[1:-1,1:])*(v[:-1,1:-1]+v[1:,1:-1]))
    # Boundary conditions
    if (BC["type"]["Left"] == "MovingWall"):
        dUVdx[0,:] = 0.5*((u[0,:-1]+u[0,1:])*(BC["v"]["Left"]))
    if (BC["type"]["Right"] == "MovingWall"):
        dUVdx[-1,:] = 0.5*((u[-1,:-1]+u[-1,1:])*(BC["v"]["Right"]))
    return dUVdx

def D2VDY2(v,mesh):
    return (v[:,1:]-v[:,:-1])/mesh.dy

def D2VDX2(v,mesh,BC):
    d2vdx2 = np.zeros((mesh.Nx+1,mesh.Ny-1))
    d2vdx2[1:-1,:] = (v[1:,1:-1]-v[:-1,1:-1])/mesh.dx
    # Boundary conditions
    if (BC["type"]["Left"] == "MovingWall"):
        d2vdx2[0,:] = 2*(v[0,1:-1]-BC["v"]["Left"])/mesh.dx
    elif (BC["type"]["Left"] == "Wall"):
        d2vdx2[0,:] = 2*v[0,1:-1]/mesh.dx
    if (BC["type"]["Right"] == "MovingWall"):
        d2vdx2[-1,:] = 2*(BC["v"]["Right"]-v[-1,1:-1])/mesh.dx
    elif (BC["type"]["Right"] == "Wall"):
        d2vdx2[-1,:] = -2*v[-1,1:-1]/mesh.dx      
    return d2vdx2

def DUDX(u,mesh):
    return u

def DVDY(v,mesh):
    return v

def STRESS(u,v,mesh,BC):
    d2udx2 = D2UDX2(u,mesh)
    d2udy2 = D2UDY2(u,mesh,BC)
    d2vdy2 = D2VDY2(v,mesh)
    d2vdx2 = D2VDX2(v,mesh,BC)
    return 1./mesh.dx*(d2udx2[1:,:]-d2udx2[:-1,:])+1./mesh.dy*(d2udy2[:,1:]-d2udy2[:,:-1]),1./mesh.dy*(d2vdy2[:,1:]-d2vdy2[:,:-1])+1./mesh.dx*(d2vdx2[1:,:]-d2vdx2[:-1,:])

def TRANS(u,v,mesh,BC):
    duudx = DUUDX(u,mesh)
    duvdy = DUVDY(u,v,mesh,BC)
    dvvdy = DVVDY(v,mesh)
    duvdx = DUVDX(u,v,mesh,BC)
    return -1./mesh.dx*(duudx[1:,:]-duudx[:-1,:])-1./mesh.dy*(duvdy[:,1:]-duvdy[:,:-1]), -1./mesh.dy*(dvvdy[:,1:]-dvvdy[:,:-1])-1./mesh.dx*(duvdx[1:,:]-duvdx[:-1,:])

def DIV(u,v,mesh):
    dudx = DUDX(u,mesh)
    dvdy = DVDY(v,mesh)
    return 1./mesh.dx*(dudx[1:,:]-dudx[:-1,:]) + 1./mesh.dy*(dvdy[:,1:]-dvdy[:,:-1])


def GRAD(p,mesh,BC):
    gradx = np.zeros_like(mesh.Xex)
    grady = np.zeros_like(mesh.Xey)
    gradx[1:-1,:] = (p[1:,:]-p[:-1])/mesh.dx
    grady[:,1:-1] = (p[:,1:] - p[:,:-1])/mesh.dy
    # Boundary conditions
    #if BC["type"]["Left"] == "Outflow":
    #    gradx[0,:] = (p[1,:]-BC["p"]["Left"])/mesh.dx
    #if BC["type"]["Right"] == "Outflow":
    #    gradx[-1,:] = (BC["p"]["Right"]-p[-2,:])/mesh.dx
    # Nothing to do, Neumann in p
    return gradx, grady 
