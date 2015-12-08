#!/usr/bin/python

#####        Labrosse 2001/2003/2015
#####     Age and evolution of growth rate of the inner core



import numpy as np
import scipy.integrate as scODE
import matplotlib.pyplot as plt



### Physical parameters

rOC=3480.e3 #m
rICp = 1221.e3 #m

Qcmb0 = 10e12 # watts
h = 0 #radiactivity

#from fitting PREM
K0 =1403.e9 #Pa
rho0 = 12451 #kg/m-3
Kprim0 = 3.567 #no units

G = 6.67384e-11 # m3/kg/s2

DeltaS = 127 #J/K/kg
mP =9e-9 #K/Pa
mX = -21e3 #K
TLricp=5700 #K
X0 = 5.6e-2 

gamma=1.5

Cp=750 #J/K/kg
beta=0.83


dTsdTad=1.65
k=200.#W/K/m
kappa=k/rho0/Cp
print "\kappa = ", kappa

gprim=4./rICp

### Lenghts and constants
Lrho =np.sqrt( 3.*K0/(2.*np.pi *G *rho0**2.) )
Arho =(5.*Kprim0-13.)/10.

    

    

### Numerical parameters

Nc=100 #number of steps in c (radius IC)




### Useful functions

def calcDensity(r):
    """ Calculate the density at the radius r """
    return rho0*(1.-r**2./Lrho**2.-Arho*r**4./Lrho**4.)

def calcGravity(r):
    """ Calculate the gravity at the radius r """
    return 4.*np.pi/3.*G*rho0*r*(1-3./5.*r**2./Lrho**2.-3.*Arho/7.*r**4./Lrho**4.)
    
def calcFC(x,d):
    """ Calculate the approximate value of function fC(x,d) (as defined in Labrosse 2015 eq. A.1.) """
    return x**3.* ( 1 - 3./5. * (d+1) * x**2. - 3./14. *(d+1) * (2*Arho-d) * x**4. )

def calcFxi(x,r_):
    """ return the function f_\chi, equation A.15 in Labrosse 2015"""
    return x**3. * ( -r_**2/3./Lrho**2. + 1./5. * (1+r_**2/Lrho**2.)*x**2. - 13./70. *x**4.)


TL0=TLricp + K0*mP*rICp**2./Lrho**2. + mX*X0*rICp**3./(Lrho**3.*fcOC)


def calcTempFusion(r):
    """ Calculate the melting temperature at the IC boundary at rIC """
    return TL0 -K0*mP*r**2./Lrho**2. + mX*X0*r**3./(Lrho**3.*calcFC(rOC/Lrho,0.))


def calcDTLDrIC(r):
    """ Calculate the diff of melting temperature at the IC boundary at rIC """
    return -K0*mP*2.*r/Lrho**2. + 3.*mX*X0*r**2./(Lrho**3.*calcFC(rOC/Lrho,0.))

def calcA7intFunction(r):
    """ Voir equation A.7 in Labrosse 2015, this is the function inside the integrale """
    return r**2.*calcDensity(r)**(gamma+1.)

def calcA13intFunction(x,ric_):
    """ Voir equation A.13 in Labrosse 2015, this is the function inside the integrale """
    A=x**2.*(1-x**2.-Arho*x**4.)*(x**2-ric_**2./Lrho**2)*(1.-3./10.*(x**2.+ric_**2./Lrho**2.))
    return A


##  Functions for the P (PL, PC, PX)

def calcPC(r_):
    """ from equation A7 (Labrosse 2015) """
    result, err=scODE.quad(calcA7intFunction,r_,rOC)
    Pc=-4.*np.pi * Cp/ calcDensity(r_)**gamma\
      *( calcDTLDrIC(r_) + 2.*gamma * rho0 *calcTempFusion(r_) *r_ / (calcDensity(r_)*Lrho**2.)*(1+2.*Arho*r_**2./Lrho**2.)  ) *result
    return Pc

def calcPC2(r):
    """ from equation A.8 (Labrosse 2015) """
    Pc2=-4.*np.pi/3.*rho0*Cp*Lrho**3. \
      *( 1-r**2./Lrho**2-Arho*r**4./Lrho**4. )**(-gamma)\
      *(calcDTLDrIC(r)+2.*gamma*calcTempFusion(r)*r/Lrho**2.*(1+2.*Arho*r**2./Lrho**2)/ ( 1-r**2./Lrho**2.-Arho*r**4./Lrho**4. ) ) \
      *(calcFC(rOC/Lrho,gamma)-calcFC(r/Lrho,gamma))
    return Pc2

def calcPL(r):
    """ from equation A.5 (Labrosse 2015) """
    return 4.*np.pi*r**2.*calcTempFusion(r)*calcDensity(r)*DeltaS


def calcPX(r):
    """ from equation A.13 (Labrosse 2015)"""
    result, err=scODE.quad(calcA13intFunction,r/Lrho,rOC/Lrho,args=r)
    return 8.*np.pi**2.*X0*G*rho0**2.*beta*r**2.*Lrho**2/calcFC(rOC/Lrho,0)*result

def calcPX2(r):
    """ from equation A.14 (Labrosse 2015) """
    return 8*np.pi**2 *X0 *G *rho**2*beta *r**2.*Lrho**2./calcFC(rOC/Lrho,0)*(calcFxi(rOC/Lrho,r)-calcFxi(r/Lrho,r))


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



r=rICp  ### Calcul pour le dernier pas de temps


### Latent heat

mathcalL_approx=4.*np.pi/3. *rho0 * TL0*DeltaS * r**3. \
  * (  1-3./5.*( 1+K0/TL0*mP )*r**2./Lrho**2.*X0/(2*calcFC(rOC/Lrho,0.)*TL0)*mX*r**3./Lrho**3. )

mathcalL,err=scODE.quad(calcPL,0,rICp)

print "latent heat", mathcalL,mathcalL_approx

### Secular cooling

mathcalC_approx=4.*np.pi/3.*rho0*Cp*Lrho*r**2*calcFC(rOC/Lrho,gamma)\
  *( mP*K0-gamma*TL0-mX*X0/calcFC(rOC/Lrho,0.)*r/Lrho )

mathcalC,err=scODE.quad(calcPC,0,rICp)

print "secular cooling", mathcalC, mathcalC_approx

### Compositional energy

mathcalX,err=scODE.quad(calcPX,0,rICp)
mathcalX2,err=scODE.quad(calcPX2,0,rICp)
print "compositional energy", mathcalX, mathcalX2

print "Total energy", mathcalX+ mathcalC+mathcalL

### Age IC

Aic=(mathcalL+mathcalX+mathcalC)/Qcmb0
print 'Qcmb = ', Qcmb0, ' Watts'
print Aic/(np.pi*1e7)/1e9, ' Gyrs'
print Qcmb0/(calcPL(r)+calcPC(r)+calcPX(r)), ' m/s'


plt.plot(np.linspace(7e12,15e12,20),(mathcalL+mathcalX+mathcalC)/np.linspace(7e12,15e12,20)/(np.pi*1e7)/1e9)


### Graphs

Qcmb=Qcmb0+np.arange(-3,7,2)*1e12
print Qcmb

f, axarr = plt.subplots(2,1)
f2, axarr2 = plt.subplots(2,1)
f3, axarr3 = plt.subplots(2,1)

c=np.linspace(0.1*rICp,rICp,200)

for Q in Qcmb:
    
    
    dcdt=np.zeros(200)
    t=np.zeros(200)
    Tic=np.zeros(200)
    S=np.zeros(200)
    
    for i in range(0,200):
        dcdt[i] = Q/(calcPL(c[i])+calcPC(c[i])+calcPX(c[i]))
        t[i] = (calcmathcalL(c[i])+calcmathcalX(c[i])+calcmathcalC(c[i]))/Q/(np.pi*1.e7*1.e6)
        Tic[i] = 3. * kappa/ (dTsdTad-1) * ( calcPL(c[i])+calcPC(c[i])+calcPX(c[i]) )/ Q /c[i]
        S[i] = 3* kappa * rho0 * gprim * gamma * TL0/K0 *( 1./Tic[i]-1)

    tauIC=(mathcalL+mathcalX+mathcalC)/Q/(np.pi*1.e7*1.e6)
   
    axarr[0].plot(c/1.e3,dcdt)
    axarr[1].plot(t,c/1.e3)
    axarr[1].plot(t,(rICp)*(t/tauIC)**0.4/1.e3)
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
        Tic[i] = 3. * kappa/ (dTsdTad-1) * ( calcPL(c[i])+calcPC2(c[i])+calcPX2(c[i]) )/ (mathcalX+ mathcalC+mathcalL) /c[i] *tauIC
        S[i] = 3* kappa * rho0 * gprim * gamma * TLricp/K0 *( 1./Tic[i]-1)


    axarr3[0].plot(c/1.e3,Tic)
    axarr3[0].plot(c/1.e3,np.ones(200))
    axarr3[1].plot(c/1.e3,S*np.pi*1e7*1e9)
    axarr3[1].plot(c/1.e3,np.zeros(200))

axarr3[0].set_title('Tic et S as fn of radius for diff age of IC')
print TL0, calcTempFusion(rICp)



Pl=np.zeros(200)
Pc=np.zeros(200)
Px=np.zeros(200)

for i in range(0,200):
    Pl[i] = calcmathcalL(c[i])
    Px[i] = calcmathcalX(c[i])
    Pc[i] = calcmathcalC(c[i])

t=(Pl+Px+Pc)/Qcmb0/(np.pi*1.e7*1.e6)

f4, axarr4 = plt.subplots(1,2)
axarr4[0].plot(c/1.e3, Pl, c/1.e3, Px,c/1.e3, Pc,c/1.e3, Pl+Px+Pc)
axarr4[1].plot(t, Pl, t, Px,t, Pc,t, Pl+Px+Pc)
axarr4[0].set_title('energy as fn of radius (km)')
axarr4[1].set_title('energy as fn of time (Ma)')

plt.show()
