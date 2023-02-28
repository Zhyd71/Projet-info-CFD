# This program is developped for academic purposes 
# GNU General Public License
# copyright: jorge.brandle@coria.fr

import numpy as np
def init_c(f,mesh):
    c = np.zeros_like(mesh.triang.triangles[:,0])*1.0
    for i, tri in enumerate(mesh.triang_fine.triangles):
        if (f(mesh.xcellcen_fine[i],mesh.ycellcen_fine[i])):
            c[mesh.father[i]] +=  1.0/4**mesh.Nrefin
    return c

def bande(x,y):
    return abs(x-50)<30

def bande_h(x,y):
    return abs(y-50)<30

def cercle(x,y):
    return np.sqrt((x-0.)**2+(y-0.0)**2)<15

def zalesak(x,y):
    return cercle(x,y)*(1.0-(x<0+2.5)*(x>0-2.5)*(y<0+10))

def one(x,y):
    return 1.0
