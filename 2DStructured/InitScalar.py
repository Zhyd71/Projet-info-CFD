# This program is developped for academic purposes 
# GNU General Public License
# copyright: jorge.brandle@coria.fr

import numpy as np
def init_c(f,mesh):
    ### Simple method not very accurate
    #c = f(mesh.Xn,mesh.Yn)

    c = (f(mesh.Xn+0.25*mesh.dx,mesh.Yn+0.25*mesh.dy)+
         f(mesh.Xn-0.25*mesh.dx,mesh.Yn+0.25*mesh.dy)+
         f(mesh.Xn-0.25*mesh.dx,mesh.Yn-0.25*mesh.dy)+
         f(mesh.Xn+0.25*mesh.dx,mesh.Yn-0.25*mesh.dy))/4.0

    ### Too expensive method
    #c = np.zeros_like(mesh.Xn)*1.0
    #Npoints = 5
    #for i in xrange(mesh.Nx):
    #    for j in xrange(mesh.Ny):
    #        for ii in xrange(Npoints):
    #            x = mesh.Xn[i,j] + mesh.dx*((1.*ii)/Npoints-0.5)
    #            for jj in xrange(Npoints):
    #                y = mesh.Yn[i,j] + mesh.dy*((1.*jj)/Npoints-0.5)
    #                if (f(x,y)):
    #                   c[i,j] = c[i,j] + 1.0/(Npoints*Npoints)
    return c

def bande(x,y):
    return 1.0*(np.abs(x-50)<30)

def bande_h(x,y):
    return 1.0*(np.abs(y-50)<30)

def cercle(x,y):
    return 1.0*(np.sqrt((x-0.)**2+(y-0.0)**2)<15)

def sbande(x,y):
    return 1.0*(np.abs(x-0.)<0.15)

def sbande_h(x,y):
    return 1.0*(np.abs(y-0.)<0.15)

def scercle(x,y):
    return 1.0*(np.sqrt((x-0.)**2+(y-0.0)**2)<0.15)

def zalesak(x,y):
    return cercle(x,y)*1.0*((1.0-(x<0+2.5)*(x>0-2.5)*(y<0+10)))

def one(x,y):
    return 1.0
