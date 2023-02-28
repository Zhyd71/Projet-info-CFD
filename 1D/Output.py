# This program is developped for academic purposes 
# GNU General Public License
# copyright: jorge.brandle@coria.fr

import matplotlib.pyplot as plt

def savefigC(c,t,x,clim=(0,1)):
    plt.cla()
    plt.ylim(clim)
    plt.plot(x,c,"b-o")
    plt.title("t="+'{:07.3f}'.format(t))
    plt.xlabel("x")
    plt.ylabel("c")
    plt.savefig("Solution/c"+'{:07.3f}'.format(t)+".png")
    plt.close()

def savefigEuler(rho,rhou,rhoE,p,c,t,x):
    plt.cla()
    plt.clf()
    plt.title("t="+'{:07.3f}'.format(t))

    plt.subplot(321)             
    plt.plot(x,rho,"b-o")
    plt.ylabel("rho")

    plt.subplot(323)             
    plt.plot(x,rhou,"b-o")
    plt.ylabel("rhou")

    plt.subplot(325)             
    plt.plot(x,rhoE,"b-o")
    plt.ylabel("rhoE")

    plt.subplot(322)             
    plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
    plt.plot(x,c,"b-o")
    plt.ylabel("c")

    plt.subplot(324)             
    plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
    plt.plot(x,rhou/rho,"b-o")
    plt.ylabel("u")

    plt.subplot(326)             
    plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
    plt.plot(x,p,"b-o")
    plt.ylabel("p")

    plt.savefig("Solution/Euler"+'{:07.3f}'.format(t)+".png")
    plt.close()



