""" Code to look at thermal evolution of a core """


import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt



class Simulation(object):

    def __init__(self):
        pass

    



class Profile(object):

    def __init__(self):
        self.Rc = 3500e3 # radius core
        self.properties = Properties()


    def T_ad(self, r, Tc):
        "adiabate, anchored at the CMB"
        return Tc * np.exp((r**2-self.Rc**2)/self.properties.LT**2)

    def rho(self, r):
        return self.properties.rhoc

    def int_ad(self, r0, r1):
        def to_integrate(r):
            return self.T_ad(r,1.)*self.rho(r)*4*np.pi*r**2
        return self.properties.Cp * scipy.integrate.quad(to_integrate, r0, r1)[0]



class Properties(object):

    def __init__(self):
        self.LT = 6400e3 # typical length adiabat
        self.rhoc = 0.8e4 # average density
        self.Cp = 750.

    def T_liquidus(self, P, xi):
        pass





def evolution_before_F():
    run = Profile()
    T_init = 6000
    Qcmb = 10e12 # 10 tera watts
    print(run.int_ad(0., run.Rc))
    dTdt = Qcmb  / run.int_ad(0., run.Rc)

    print(dTdt*np.pi*1e7*1e9)

# -4*np.pi*run.Rc**2

if __name__ == "__main__":

    evolution_before_F()