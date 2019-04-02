""" Code to look at thermal evolution of a core """


import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt



class Simulation(object):

    def __init__(self):
        pass #initialisation

    def run(self):
        pass #run 

    def run_1step(self):
        pass # run 1 step (to be used in self.run())



class Profile(object):

    def __init__(self):
        self.Rc = 3500e3 # radius core
        self.Qcmb = 10e12 # 10 terawatts
        self.XOC = 3.6 #weigth percents
        self.properties = Properties()


    def T_ad(self, r, Tc):
        "adiabat temperature, anchored at the CMB"
        return Tc * np.exp((r**2-self.Rc**2)/self.properties.LT**2)

    def rho(self, r):
        return self.properties.rhoc


    def pressure(self, r):
        pass

    def int_ad(self, r0, r1):
        def to_integrate(r):
            return self.T_ad(r,1.)*self.rho(r)*4*np.pi*r**2
        return self.properties.Cp * scipy.integrate.quad(to_integrate, r0, r1)[0]

    def evolution_before_F(self):
        dTdt = self.Qcmb  / self.int_ad(0., self.Rc)
        return dTdt

class Properties(object):

    def __init__(self):
        self.LT = 6400e3 # typical length adiabat
        self.rhoc = 0.8e4 # average density
        self.Cp = 750.
        self.dTliqdP = 10e-9 #K/Pa
        self.dTliqdX = 6.3 #K/weigth percents
        self.Teu = 6000 # eutectic temperature (K) at pressure P_0
        self.Xeu = 1.5 #eutectic composition

    def T_liquidus(self, P, xi):
        P_0 = 250e9 # ? Pressure ICB?
        return self.Teu + (P-P_0)*self.dTliqdP + (xi-self.Xeu)*dTliqdX
    






# -4*np.pi*run.Rc**2

if __name__ == "__main__":

    evolution_before_F()