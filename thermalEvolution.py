#!/usr/bin/python

#####        Labrosse 2001/2003/2015
#####     Age and evolution of growth rate of the inner core



import numpy as np
import scipy.integrate as scODE
import matplotlib.pyplot as plt



### Physical parameters

R_OC = 3480.e3 #m
R_IC_P = 1221.e3 #m

QCMB = 10e12 # watts
H = 0 #radiactivity  !!! Not implemented yet

#from fitting PREM
K0 = 1403.e9 #Pa
RHO_0 = 12451 #kg/m-3
K_PRIM0 = 3.567 #no units

G = 6.67384e-11 # m3/kg/s2

DELTA_S = 127 #J/K/kg
M_P = 9e-9 #K/Pa
M_X = -21e3 #K
TL_RIC = 5700 #K
X0 = 5.6e-2 
GAMMA = 1.5
CP = 750 #J/K/kg
BETA = 0.83


DTS_DTAD = 1.65
k = 200.#W/K/m
KAPPA = k/RHO_0/CP
print "\kappa = ", KAPPA

G_PRIM = 4./R_IC_P

### Lenghts and constants
L_RHO = np.sqrt(3.*K0/(2.*np.pi *G *RHO_0**2.))
A_RHO = (5.*K_PRIM0-13.)/10.

    

    

### Numerical parameters

Nc = 100 #number of steps in c (radius IC)




### Useful functions

def calcDensity(r):
    """ Calculate the density at the radius r """
    return RHO_0*(1.-r**2./L_RHO**2.-A_RHO*r**4./L_RHO**4.)

def calcGravity(r):
    """ Calculate the gravity at the radius r """
    return 4.*np.pi/3.*G*RHO_0*r*(1-3./5.*r**2./L_RHO**2.-3.*A_RHO/7.*r**4./L_RHO**4.)
    
def calcFC(x,d):
    """ Calculate the approximate value of function fC(x,d) (as defined in Labrosse 2015 eq. A.1.) """
    return x**3.* ( 1 - 3./5. * (d+1) * x**2. - 3./14. *(d+1) * (2*A_RHO-d) * x**4. )

def calcFxi(x,r_):
    """ return the function f_\chi, equation A.15 in Labrosse 2015"""
    return x**3. * ( -r_**2/3./L_RHO**2. + 1./5. * (1+r_**2/L_RHO**2.)*x**2. - 13./70. *x**4.)




def calcTempFusion(r):
    """ Calculate the melting temperature at the IC boundary at rIC """
    return TL0 -K0*M_P*r**2./L_RHO**2. + M_X*X0*r**3./(L_RHO**3.*calcFC(R_OC/L_RHO,0.))


def calcDTLDrIC(r):
    """ Calculate the diff of melting temperature at the IC boundary at rIC """
    return -K0*M_P*2.*r/L_RHO**2. + 3.*M_X*X0*r**2./(L_RHO**3.*calcFC(R_OC/L_RHO,0.))

def calcA7intFunction(r):
    """ Voir equation A.7 in Labrosse 2015, this is the function inside the integrale """
    return r**2.*calcDensity(r)**(GAMMA+1.)

def calcA13intFunction(x,ric_):
    """ Voir equation A.13 in Labrosse 2015, this is the function inside the integrale """
    A=x**2.*(1-x**2.-A_RHO*x**4.)*(x**2-ric_**2./L_RHO**2)*(1.-3./10.*(x**2.+ric_**2./L_RHO**2.))
    return A


##  Functions for the P (PL, PC, PX)

def calcPC(r_):
    """ from equation A7 (Labrosse 2015) """
    result, err=scODE.quad(calcA7intFunction,r_,R_OC)
    Pc=-4.*np.pi * CP/ calcDensity(r_)**GAMMA\
      *( calcDTLDrIC(r_) + 2.*GAMMA * RHO_0 *calcTempFusion(r_) *r_ / (calcDensity(r_)*L_RHO**2.)*(1+2.*A_RHO*r_**2./L_RHO**2.)  ) *result
    return Pc

def calcPC2(r):
    """ from equation A.8 (Labrosse 2015) """
    Pc2=-4.*np.pi/3.*RHO_0*CP*L_RHO**3. \
      *( 1-r**2./L_RHO**2-A_RHO*r**4./L_RHO**4. )**(-GAMMA)\
      *(calcDTLDrIC(r)+2.*GAMMA*calcTempFusion(r)*r/L_RHO**2.*(1+2.*A_RHO*r**2./L_RHO**2)/ ( 1-r**2./L_RHO**2.-A_RHO*r**4./L_RHO**4. ) ) \
      *(calcFC(R_OC/L_RHO,GAMMA)-calcFC(r/L_RHO,GAMMA))
    return Pc2

def calcPL(r):
    """ from equation A.5 (Labrosse 2015) """
    return 4.*np.pi*r**2.*calcTempFusion(r)*calcDensity(r)*DELTA_S


def calcPX(r):
    """ from equation A.13 (Labrosse 2015)"""
    result, err=scODE.quad(calcA13intFunction,r/L_RHO,R_OC/L_RHO,args=r)
    return 8.*np.pi**2.*X0*G*RHO_0**2.*BETA*r**2.*L_RHO**2/calcFC(R_OC/L_RHO,0)*result

def calcPX2(r):
    """ from equation A.14 (Labrosse 2015) """
    return 8*np.pi**2 *X0 *G *rho**2*BETA *r**2.*L_RHO**2./calcFC(R_OC/L_RHO,0)*(calcFxi(R_OC/L_RHO,r)-calcFxi(r/L_RHO,r))


## Functions pour mathcal L, C, X (valeurs les plus precises)

def calcmathcalL(r):
    results,err=scODE.quad(calcPL,0,r)
    return results

def calcmathcalC(r):
    results,err=scODE.quad(calcPC,0,r)
    return results

def calcmathcalX(r):
    results,err=scODE.quad(calcPX,0,r)
    return results

#########
#########



r = R_IC_P  ### Calcul pour r=r_ic

fcOC = 
TL0 = TL_RIC + K0*M_P*R_IC_P**2./L_RHO**2. + M_X*X0*R_IC_P**3./(L_RHO**3.*fcOC)


### Latent heat

mathcalL_approx=4.*np.pi/3. *RHO_0 * TL0*DELTA_S * r**3. \
  * (  1-3./5.*( 1+K0/TL0*M_P )*r**2./L_RHO**2.*X0/(2*calcFC(R_OC/L_RHO,0.)*TL0)*M_X*r**3./L_RHO**3. )

mathcalL,err=scODE.quad(calcPL,0,R_IC_P)

print "latent heat", mathcalL,mathcalL_approx

### Secular cooling

mathcalC_approx=4.*np.pi/3.*RHO_0*CP*L_RHO*r**2*calcFC(R_OC/L_RHO,GAMMA)\
  *( M_P*K0-GAMMA*TL0-M_X*X0/calcFC(R_OC/L_RHO,0.)*r/L_RHO )

mathcalC,err=scODE.quad(calcPC,0,R_IC_P)

print "secular cooling", mathcalC, mathcalC_approx

### Compositional energy

mathcalX,err=scODE.quad(calcPX,0,R_IC_P)
mathcalX2,err=scODE.quad(calcPX2,0,R_IC_P)
print "compositional energy", mathcalX, mathcalX2

print "Total energy", mathcalX+ mathcalC+mathcalL

### Age IC

Aic=(mathcalL+mathcalX+mathcalC)/QCMB
print 'Qcmb = ', QCMB, ' Watts'
print Aic/(np.pi*1e7)/1e9, ' Gyrs'
print QCMB/(calcPL(r)+calcPC(r)+calcPX(r)), ' m/s'


plt.plot(np.linspace(7e12,15e12,20),(mathcalL+mathcalX+mathcalC)/np.linspace(7e12,15e12,20)/(np.pi*1e7)/1e9)


### Graphs

Qcmb=QCMB+np.arange(-3,7,2)*1e12
print Qcmb

f, axarr = plt.subplots(2,1)
f2, axarr2 = plt.subplots(2,1)
f3, axarr3 = plt.subplots(2,1)

c=np.linspace(0.1*R_IC_P,R_IC_P,200)

for Q in Qcmb:
    
    
    dcdt=np.zeros(200)
    t=np.zeros(200)
    Tic=np.zeros(200)
    S=np.zeros(200)
    
    for i in range(0,200):
        dcdt[i] = Q/(calcPL(c[i])+calcPC(c[i])+calcPX(c[i]))
        t[i] = (calcmathcalL(c[i])+calcmathcalX(c[i])+calcmathcalC(c[i]))/Q/(np.pi*1.e7*1.e6)
        Tic[i] = 3. * KAPPA/ (DTS_DTAD-1) * ( calcPL(c[i])+calcPC(c[i])+calcPX(c[i]) )/ Q /c[i]
        S[i] = 3* KAPPA * RHO_0 * G_PRIM * GAMMA * TL0/K0 *( 1./Tic[i]-1)

    tauIC=(mathcalL+mathcalX+mathcalC)/Q/(np.pi*1.e7*1.e6)
   
    axarr[0].plot(c/1.e3,dcdt)
    axarr[1].plot(t,c/1.e3)
    axarr[1].plot(t,(R_IC_P)*(t/tauIC)**0.4/1.e3)
    axarr2[0].plot(c/1.e3,Tic)
    axarr2[0].plot(c/1.e3,np.ones(200))
    axarr2[1].plot(c/1.e3,S*np.pi*1e7*1e9)
    axarr2[1].plot(c/1.e3,np.zeros(200))

axarr[0].set_title('dcdt as fn of c')


axarr[1].set_title('c as fn of t')
axarr2[0].set_title('Tic and S as fn of radius for different values of Q')

   
for tauIC in (np.arange(200,1800,200)*1.e6*1.e7*np.pi):

    Tic=np.zeros(200)
    S=np.zeros(200)

    for i in range(0,200):
        Tic[i] = 3. * KAPPA/ (DTS_DTAD-1) * ( calcPL(c[i])+calcPC2(c[i])+calcPX2(c[i]) )/ (mathcalX+ mathcalC+mathcalL) /c[i] *tauIC
        S[i] = 3* KAPPA * RHO_0 * G_PRIM * GAMMA * TL_RIC/K0 *( 1./Tic[i]-1)


    axarr3[0].plot(c/1.e3,Tic)
    axarr3[0].plot(c/1.e3,np.ones(200))
    axarr3[1].plot(c/1.e3,S*np.pi*1e7*1e9)
    axarr3[1].plot(c/1.e3,np.zeros(200))

axarr3[0].set_title('Tic et S as fn of radius for diff age of IC')
print TL0, calcTempFusion(R_IC_P)



Pl=np.zeros(200)
Pc=np.zeros(200)
Px=np.zeros(200)

for i in range(0,200):
    Pl[i] = calcmathcalL(c[i])
    Px[i] = calcmathcalX(c[i])
    Pc[i] = calcmathcalC(c[i])

t=(Pl+Px+Pc)/QCMB/(np.pi*1.e7*1.e6)

f4, axarr4 = plt.subplots(1,2)
axarr4[0].plot(c/1.e3, Pl, c/1.e3, Px,c/1.e3, Pc,c/1.e3, Pl+Px+Pc)
axarr4[1].plot(t, Pl, t, Px,t, Pc,t, Pl+Px+Pc)
axarr4[0].set_title('energy as fn of radius (km)')
axarr4[1].set_title('energy as fn of time (Ma)')

plt.show()
