# This program is developped for academic purposes 
# GNU General Public License
# copyright: jorge.brandle@coria.fr

import numpy as np


def InitRest(mesh):
    # Fluid at rest
    return np.zeros_like(mesh.Xex), np.zeros_like(mesh.Xey), np.zeros_like(mesh.Xn)
    
