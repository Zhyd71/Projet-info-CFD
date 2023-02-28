# This program is developped for academic purposes 
# GNU General Public License
# copyright: jorge.brandle@coria.fr

import numpy as np
from matplotlib.tri import Triangulation, UniformTriRefiner

def Compute_CellCen(triang):
    xcellcen = np.zeros(np.size(triang.triangles[:,0]))
    ycellcen = np.zeros(np.size(triang.triangles[:,0]))
    for i,tri in enumerate(triang.triangles):
        xcellcen[i] = (triang.x[tri[0]]+triang.x[tri[1]]+triang.x[tri[2]])/3.
        ycellcen[i] = (triang.y[tri[0]]+triang.y[tri[1]]+triang.y[tri[2]])/3.
    triangcellcent = Triangulation(xcellcen, ycellcen)
    return xcellcen, ycellcen, triangcellcent

def Compute_EdgeCen(triang):
    xedgecen = np.zeros(np.size(triang.edges[:,0]))
    yedgecen = np.zeros(np.size(triang.edges[:,0]))
    for i,edge in enumerate(triang.edges):
        xedgecen[i] = (triang.x[edge[0]] + triang.x[edge[1]])/2
        yedgecen[i] = (triang.y[edge[0]] + triang.y[edge[1]])/2
    return xedgecen, yedgecen

def PointInTriangle(xp,yp,x1,y1,x2,y2,x3,y3):
    def sign(x1,y1,x2,y2,x3,y3):
        return (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3)
    b1 = sign(xp,yp,x1,y1,x2,y2) < 0.0
    b2 = sign(xp,yp,x2,y2,x3,y3) < 0.0
    b3 = sign(xp,yp,x3,y3,x1,y1) < 0.0

    return ((b1 == b2) and (b2 == b3))

def ComputeArea(triang):
    Area = np.zeros(np.size(triang.triangles[:,0]))
    for i,tri in enumerate(triang.triangles):
        Area[i] = 0.5 * ((triang.x[tri[0]]-triang.x[tri[2]])*(triang.y[tri[1]]-triang.y[tri[0]])-
                         (triang.x[tri[0]]-triang.x[tri[1]])*(triang.y[tri[2]]-triang.y[tri[0]]))
    return Area

def ComputeFine(triang,Nrefin):
    refiner = UniformTriRefiner(triang)
    triang_fine, indexes = refiner.refine_triangulation(True,Nrefin)
    xcellcen_fine, ycellcen_fine, triangcellcent_fine = Compute_CellCen(triang_fine)
    father = np.zeros_like(xcellcen_fine,dtype=np.int16)
    for i in range(len(father)):
        father[i] = int(i/4**Nrefin)
    return xcellcen_fine, ycellcen_fine, triang_fine, father, triangcellcent_fine



# Compute normals
def ComputeNormals(triang):
    nx = np.zeros(np.size(triang.edges[:,0]))
    ny = np.zeros(np.size(triang.edges[:,0]))
    dS = np.zeros(np.size(triang.edges[:,0]))
    for i,edge in enumerate(triang.edges):
        nx[i] =   triang.y[edge[0]] - triang.y[edge[1]]
        ny[i] = -(triang.x[edge[0]] - triang.x[edge[1]])
        dS[i] = np.sqrt(nx[i]**2+ny[i]**2)
    return nx,ny,dS
    nx[i], ny[i] = nx[i]/dS[i], ny[i]/dS[i]

def FindDirections(triang):
    Ltri = -np.ones(np.size(triang.edges[:,0]),dtype=np.int16)
    Rtri = -np.ones(np.size(triang.edges[:,0]),dtype=np.int16)
    indEdge = np.zeros_like(triang.neighbors)
    direc = -1.0*np.ones_like(triang.neighbors)
#for j,tri in enumerate(triang.triangles):
#    print  j,"/",np.size(triang.triangles[:,0])	
#    for k,tri2 in enumerate(triang.neighbors[j]):
#        my_edge = np.zeros(2,dtype=np.int16)
#        indedge = 0
#        for ind in tri:
#            if ind in triang.triangles[tri2]:
#                my_edge[indedge] = ind
#                indedge += 1
#        my_edge = np.sort(my_edge)[::-1]
#        # print my_edge
#        for i,edge in enumerate(triang.edges):
#             # print my_edge, edge
#             if ((edge[0] == my_edge[0]) and (edge[1] == my_edge[1])):
#                 pass
    for i,edge in enumerate(triang.edges):
    #print i,"/",np.size(triang.edges[:,0])
        for j,tri in enumerate(triang.triangles):
            if ((edge[0] in tri) and (edge[1] in tri)):
               P1 = np.nonzero(tri==edge[0])[0][0]
               P2 = np.nonzero(tri== edge[1])[0][0]
               if ( ((P1==0)and(P2==1)) or ((P1==1)and(P2==2)) or ((P1==2)and(P2==0))):
                   Ltri[i] = j
               else:
                   Rtri[i] = j          
               for k,tri2 in enumerate(triang.neighbors[j]):
                   if ((edge[0] in triang.triangles[tri2]) and (edge[1] in triang.triangles[tri2])):
                      indEdge[j,k] = i
                      if ( ((P1==0)and(P2==1)) or ((P1==1)and(P2==2)) or ((P1==2)and(P2==0))):
                          direc[j,k] = 1.0
    return Ltri, Rtri, indEdge, direc
