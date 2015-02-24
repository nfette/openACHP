# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 22:25:20 2015

@author: nfette
"""

from CoolProp.CoolProp import PropsSI
from hw2_1 import CelsiusToKelvin as C2K

MW_LiBr = 86.85 # kg/kmol
MW_H2O = 18.015 # kg/kmol

def massfraction_LiBrH2O(x):
    """input: mole fraction, x, of LiBr"""
    return x * MW_LiBr / (x * MW_LiBr + (1 - x) * MW_H2O)
    
def molefraction_LiBrH2O(w):
    """input: mass fraction, w, of LiBr"""
    return (w / MW_LiBr) / (w/MW_LiBr + (1 - w) / MW_H2O)
    
def P_LiBrH2O(T,x):
    """Return pressure above a lithium bromide-water mixture
    at the given temperature and composition using the formulation
    presented by
    Patek and Klomfar, Int. J. Refrig., Vol 29, pp 566-578 (2006)
    
    Notes: "above" the mixture: is completely water vapor. So there are only
    two relevant properties to find equilibrium vapor pressure (?).

    Units: T [K]
           x = mass fraction LiBr
           P [bar]
    """
    a = [-2.41303e2, 1.91750e7, -1.75521e8, 3.25432e7, 3.92571e2, -2.12626e3, 1.85127e8, 1.91216e3]
    m = [3,4,4,8,1,1,4,6]
    n = [0,5,6,3,0,2,6,0]
    t = [0,0,0,0,1,1,1,1]
    T_c = 647.096 # K
    T_0 = 221 # K
    h_c = 37548.5 # [J/gmol] "Enthalpy at the critical point for pure water"
    #TK = convertTemp(T$,K,T)	"convert to K"
    TK = T
    molef = x
    #if (UnitSystem('mass')=1) then x=molefraction_LiBrH2O(x)
    x_N = molefraction_LiBrH2O(x)
    s = 0
    for i in range(8):
        s=s+a[i]* x_N**m[i] * abs(0.4-x_N)**n[i] * (TK/T_c)**t[i]
    Theta=TK-s
    #Ts=convertTemp(K,T$,Theta)
    Ts = Theta
    #P$=UnitSystem$('Pressure')
    #P_LiBrH2O=pressure(Water,T=Ts,x=0)
    Q = 0.0
    #print("Ts, Q = {}, {}".format(Ts, Q))
    pressurePa = PropsSI('P','T',Ts,'Q',Q,'water') # Pa
    #print("pressurePa = {}".format(pressurePa)) # ()
    pressureBar = pressurePa * 1e-5
    return pressureBar
    
def T_LiBrH2O(P,x):
    """T_LiBrH2O returns the temperature of a lithium bromide-water mixture at
    the given the pressure and composition using the formulation presented by
    Patek and Klomfar, Int. J. Refrig., Vol 29, pp 566-578 (2006)
    
    Notes: "above" the mixture: is completely water vapor. So there are only
    two relevant properties to find equilibrium vapor pressure (?).

    Units: T [K]
           x = mass fraction LiBr
           P [bar]
    """
    # Just call a solver on the previously defined pressure function
    
    
if __name__ == "__main__":
    # Unit testing
    # Pressure (cf table 6.1, or rather LiBrSS7B.EES (slight difference))
    #for TC,x in ((32.7,0.567),(63.3,0.567),(89.9,0.625),(53.3,0.625),\
    #    (44.7,62.5),(77.0,0.0),(40.2,0.0),(1.3,0.0)):
    #    P = P_LiBrH2O(C2K(TC),x)
    #    print("(T,x) = ({},{}) -> P = {}".format(TC,x,P))
    for TC,x in ((50.0,0.5),):
        P = P_LiBrH2O(C2K(TC),x)
        print("(T,x) = ({},{}) -> P = {}".format(TC,x,P))
    