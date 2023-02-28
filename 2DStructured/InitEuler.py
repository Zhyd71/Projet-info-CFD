# This program is developped for academic purposes 
# GNU General Public License
# copyright: jorge.brandle@coria.fr

import numpy as np
from InitScalar import *

def CircleShock(mesh,gamma):
    U = {"rho": np.ones_like(mesh.Xn), "rhou": np.zeros_like(mesh.Xn), "rhov": np.zeros_like(mesh.Xn), "rhoE":  np.ones_like(mesh.Xn)}


    rho_l = 1.
    rho_r = 0.125
    p_l = 1.0
    p_r = 0.1
    u_l = 0.0
    u_r = 0.0
    v_l = 0.0
    v_r = 0.0
    rhou_l = rho_l*u_l
    rhou_r = rho_r*u_r
    rhov_l = rho_l*v_l
    rhov_r = rho_r*v_r
    rhoE_l = p_l/(gamma-1) + 0.5*rho_l*(u_l**2+v_l**2)
    rhoE_r = p_r/(gamma-1) + 0.5*rho_r*(u_l**2+v_l**2)

    U["rho"] = (rho_l-rho_r)*init_c(scercle, mesh)+rho_r
    U["rhou"] = (rhou_l-rhou_r)*init_c(scercle, mesh)+rhou_r
    U["rhov"] = (rhov_l-rhov_r)*init_c(scercle, mesh)+rhov_r
    U["rhoE"] = (rhoE_l-rhoE_r)*init_c(scercle, mesh)+rhoE_r

    return U

def BandeShock(mesh,gamma):
    U = {"rho": np.ones_like(mesh.Xn), "rhou": np.zeros_like(mesh.Xn), "rhov": np.zeros_like(mesh.Xn), "rhoE":  np.ones_like(mesh.Xn)}
    rho_l = 1.
    rho_r = 0.125
    p_l = 1.0
    p_r = 0.1
    u_l = 0.0
    u_r = 0.0
    v_l = 0.0
    v_r = 0.0
    rhou_l = rho_l*u_l
    rhou_r = rho_r*u_r
    rhov_l = rho_l*v_l
    rhov_r = rho_r*v_r
    rhoE_l = p_l/(gamma-1) + 0.5*rho_l*(u_l**2+v_l**2)
    rhoE_r = p_r/(gamma-1) + 0.5*rho_r*(u_l**2+v_l**2)

    U["rho"] = (rho_l-rho_r)*init_c(sbande, mesh)+rho_r
    U["rhou"] = (rhou_l-rhou_r)*init_c(sbande, mesh)+rhou_r
    U["rhov"] = (rhov_l-rhov_r)*init_c(sbande, mesh)+rhov_r
    U["rhoE"] = (rhoE_l-rhoE_r)*init_c(sbande, mesh)+rhoE_r

    return U

def Cavity(mesh):
    U = {"rho": np.ones_like(mesh.Xn), "rhou": np.zeros_like(mesh.Xn), "rhov": np.zeros_like(mesh.Xn), "rhoE":  np.ones_like(mesh.Xn)}
    Ul = {}
    Ur = {}
    Ud = {}
    Uu = {}
    for key in U:
        Ul[key] = U[key][mesh.Nx-1,:]
        Ur[key] = U[key][0,:]
        Ud[key] = U[key][:,mesh.Ny-1]
        Uu[key] = U[key][:,0]
    Uu["rhou"] = 0.01*np.ones_like(Uu["rhou"])
    Ud["rhou"] = 0.01*np.ones_like(Ud["rhou"])
    Ul["rhou"] = 0.01*np.ones_like(Ul["rhou"])
    Ur["rhou"] = 0.01*np.ones_like(Ur["rhou"])
    return U, Ul, Ur, Ud, Uu
    
