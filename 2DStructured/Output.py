# This program is developped for academic purposes 
# GNU General Public License
# copyright: jorge.brandle@coria.fr

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
def savefigC(c,t,mesh):
    fig, ax = plt.subplots()
    cmap = cm.get_cmap(name='bwr', lut=None)
#levels=levels,cmap=cmap,extend='both'
    pcol = plt.pcolormesh(mesh.Xn,mesh.Yn,c,cmap=cmap)
    cbar = fig.colorbar(pcol)
    ax.set_aspect('equal')
    plt.savefig("Solution/T"+'{:07.3f}'.format(t)+".png")
    plt.close()

def savefigEuler(U,t,mesh):
    fig, ax = plt.subplots()
    ax1 = plt.subplot(211)
    cmap = cm.get_cmap(name='bwr', lut=None)
#levels=levels,cmap=cmap,extend='both'
    pcol = plt.pcolormesh(mesh.Xn,mesh.Yn,U["rho"],cmap=cmap)
    plt.quiver(mesh.Xn[::4,::4],mesh.Yn[::4,::4],(U["rhou"]/U["rho"])[::4,::4],(U["rhov"]/U["rho"])[::4,::4])
    cbar = fig.colorbar(pcol)
    ax1.set_aspect('equal')

    ax2 = plt.subplot(212)
    cmap = cm.get_cmap(name='bwr', lut=None)
    pcol = plt.pcolormesh(mesh.Xn,mesh.Yn,U["rhoE"],cmap=cmap)
    cbar = fig.colorbar(pcol)
    ax2.set_aspect('equal')
    plt.savefig("Solution/Euler"+'{:07.3f}'.format(t)+".png")
    plt.close()

def savefigNS(u,v,p,t,mesh):
    fig, ax = plt.subplots(figsize=(11,7), dpi=100)
    ax.set_title("t="+'{:07.3f}'.format(t))
    cmap = cm.get_cmap(name='bwr', lut=None)
    pcol = plt.pcolormesh(mesh.Yn,mesh.Xn,p,cmap=cmap)
    cbar = fig.colorbar(pcol)
    un = (u[:-1,:]+u[1:,:])/2
    vn = (v[:,:-1]+v[:,1:])/2
    stepx = max(1,int(mesh.Nx/15)) # Not show all the arrows
    stepy = max(1,int(mesh.Ny/15))
    if (np.any(np.abs(un)>0) or np.any(np.abs(vn))>0):
        plt.quiver(mesh.Yn[::stepx,::stepy],mesh.Xn[::stepx,::stepy],un[::stepx,::stepy],vn[::stepx,::stepy],units='width',scale=10)
        # plt.quiver(mesh.Xn[::4,::4],mesh.Yn[::4,::4],un[::4,::4],vn[::4,::4])
    ax.set_aspect('equal')
    plt.savefig("Solution/NS"+'{:07.3f}'.format(t)+".png")
    plt.close()

    
     
