# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 22:25:20 2015

The purpose of this file is to check the compatibility of the relations from
Patel and Klomfar, Int J Refrig, Vol 20, pp 566-578 (2006) with the CoolProp
functions for water, and compare to a known set of code using these relations.

@author: nfette
"""

from CoolProp.CoolProp import PropsSI
from hw2_1 import CelsiusToKelvin as C2K, KelvinToCelsius as K2C
from scipy.optimize import minimize
import numpy as np

MW_LiBr = 86.85 # kg/kmol
MW_H2O = 18.015 # kg/kmol

def massfraction(x):
    """input: mole fraction, x, of LiBr"""
    return x * MW_LiBr / (x * MW_LiBr + (1 - x) * MW_H2O)
    
def molefraction(w):
    """input: mass fraction, w, of LiBr"""
    return (w / MW_LiBr) / (w/MW_LiBr + (1 - w) / MW_H2O)
    
def pressure(T,x):
    """Return pressure above a lithium bromide-water mixture
    at the given temperature and composition using the formulation
    presented by
    Patek and Klomfar, Int. J. Refrig., Vol 29, pp 566-578 (2006)
    
    Notes: "above" the mixture: is completely water vapor. So there are only
    two relevant properties to find equilibrium vapor pressure (?).

    Units: T [K]
           x = mass fraction LiBr
           P [bar]
           
    Based on Table 4 and Equation (1) in reference.
    """
    
    a = [-2.41303e2, 1.91750e7, -1.75521e8, 3.25432e7, 3.92571e2, -2.12626e3, 1.85127e8, 1.91216e3]
    m = [3,4,4,8,1,1,4,6]
    n = [0,5,6,3,0,2,6,0]
    t = [0,0,0,0,1,1,1,1]
    T_c = 647.096 # K
    h_c = 37548.5 # [J/gmol] "Enthalpy at the critical point for pure water"
    #TK = convertTemp(T$,K,T)	"convert to K"
    TK = T
    molef = x
    #if (UnitSystem('mass')=1) then x=molefraction_LiBrH2O(x)
    x_N = molefraction(x)
    
    s = 0
    for i in range(8):
        s=s+a[i]* x_N**m[i] * abs(0.4-x_N)**n[i] * (TK/T_c)**t[i]
    Theta=TK-s
    #Ts=convertTemp(K,T$,Theta)
    Ts = Theta
    Q = 0.0
    #print("Ts, Q = {}, {}".format(Ts, Q))
    pressurePa = PropsSI('P','T',Ts,'Q',Q,'water') # Pa
    #print("pressurePa = {}".format(pressurePa)) # ()
    pressureBar = pressurePa * 1e-5
    return pressureBar
    
def objective_T(T,*Px):
    #print("T, Px = {}, {}".format(T,Px))
    P,x = Px
    return (P - pressure(T,x)) ** 2
    
def temperature(P,x):
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
    # Trick is to get the guess temperature on the high side.
    # Does this constraint do anything? Maybe not.
    cons = ({"type": "ineq",
             "fun": lambda T: np.array(pressure(T,x) - P),
            },)
    soln = minimize(objective_T, (647.), constraints=cons, args=(P,x,))
    return soln.x[0]

def objective_x(x,*TP):
    #print("T, Px = {}, {}".format(T,Px))
    T,P = TP
    return (P - pressure(T,x)) ** 2
    
def massFraction(T,P):
    """Returns the composition of a lithium bromide-water mixture at
    the given the temprature and pressure using the formulation presented by
    Patek and Klomfar, Int. J. of Refrigeration, Vol 29, pp. 566-578, (2006)

    Notes: "above" the mixture: is completely water vapor. So there are only
    two relevant properties to find equilibrium vapor pressure (?).

    Units: T [K]
           x = mass fraction LiBr
           P [bar]    
    """
    # Just call a solver on the previously defined pressure function
    # Trick required? Guess mass fraction not too high nor low.
    # Does this constraint do anything? Maybe not.
    cons = ({"type": "ineq",
             "fun": lambda x: np.array(pressure(T,x) - P),
            },)
    soln = minimize(objective_x, (0.5), constraints=cons, args=(T,P))
    return soln.x[0]

def massSpecificEnthalpy(T,x):
    """Inputs:  T = Temperature / [Kelvin]
         x = mass fraction LiBr
Outputs: h = mass specific enthalpy / [J/kg]

Based on table 7 and equation (4) in reference.
"""
    
    a=[2.27431,-7.99511, 385.239,-16394,-422.562,0.113314,-8.33474,-17383.3,\
    6.49763,3245.52,-13464.3,39932.2,-258877,-0.00193046,2.80616,-40.4479,\
    145.342,-2.74873,-449.743,-12.1794,-0.00583739,0.233910,0.341888,8.85259,\
    -17.8731,0.0735179,-0.000179430,0.00184261,-0.00624282,0.00684765]
    m=[1,1,2,3,6,1,3,5,4,5,5,6,6,1,2,2,2,5,6,7,1,1,2,2,2,3,1,1,1,1]
    n=[0,1,6,6,2,0,0,4,0,4,5,5,6,0,3,5,7,0,3,1,0,4,2,6,7,0,0,1,2,3]
    t=[0,0,0,0,0,1,1,1,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5]
    T_crit = PropsSI('water','Tcrit') # [K]
    P_crit = PropsSI('water','pcrit')
    h_crit = PropsSI('H','T',T_crit,'P',P_crit,'water') # J/kg
    h_c = h_crit * MW_H2O
    T_c = T_crit # [K]
    T_0 = 221. # [K] "is a nonlinear parameter of the equations"

    TK = T
    x_N = molefraction(x)
    # print(x_N)
    s = 0
    for i in range(len(a)):
        s = s + a[i] * x_N**m[i] * abs(0.4-x_N)**n[i] * (T_c/(TK-T_0))**t[i]
    Qu = 0.0
    h_w_mass = PropsSI('H','T',TK,'Q',Qu,'water')
    #print("Sat liquid water enthalpy {} J/kg".format(h_w_mass))
    h_w_molar = h_w_mass * MW_H2O
    h_w_molar = h_w_molar
    h_molar = (1 - x_N) * h_w_molar + h_c * s
    MW = x_N * MW_LiBr + (1 - x_N) * MW_H2O
    #print("MW = {}".format(MW))
    result = h_molar / MW
    return result

def massSpecificEntropy(T,x):
    """Inputs:  T = Temperature / [Kelvin]
         x = mass fraction LiBr
Outputs: s = mass specific entropy / [J/kg-K]

Based on table 8 and equation (5) in reference.
"""
    
    a=[1.53091,-4.52564, 698.302,-21666.4,-1475.33,0.0847012,-6.59523,
       -29533.1,0.00956314,-0.188679,9.31752,5.78104,13893.1,-17176.2,
       415.108,-55564.7,-0.00423409,30.5242,-1.67620,14.8283,0.00303055,
       -0.0401810,0.149252,2.59240,-0.177421,-0.0000699650,0.000605007,
       -0.00165228,0.00122966]
    m = [1,1,2,3,6,1,3,5,1,2,2,4,5,5,6,6,1,3,5,7,1,1,1,2,3,1,1,1,1]
    n = [0,1,6,6,2,0,0,4,0,0,4,0,4,5,2,5,0,4,0,1,0,2,4,7,1,0,1,2,3]
    t = [0,0,0,0,0,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,4,4,4,4,4,5,5,5,5]
    T_c = 647.096 # [K]
    T_0 = 221 # [K] "is a nonlinear parameter of the equations"
    s_c = 79.3933 # [J/gmol-K]
    T_crit = PropsSI('water','Tcrit') # [K]
    P_crit = PropsSI('water','pcrit')
    s_crit = PropsSI('S','T',T_crit,'P',P_crit,'water') # J/kg-K
    s_crit_molar = s_crit * MW_H2O
    print("By the way, pure water @ critical point has entropy {} J/kg-K"
        .format(s_crit))
    print("By the way, pure water @ critical point has entropy {} J/kmol-K"
        .format(s_crit_molar))
        
    TK = T
    x_N = molefraction(x)
    s = 0
    for i in range(len(a)):
         s = s + a[i] * x_N ** m[i] * abs(0.4 - x_N) ** n[i] * (T_c / (TK - T_0)) ** t[i]
    #s_w=entropy(Water,T=T,x=0)
    Qu_water = 0.0
    s_w_mass = PropsSI('S','T',TK,'Q',Qu_water,'water') # J/kg-K
    s_w_molar = s_w_mass * MW_H2O # J/kmol-K
    s_molar = (1 - x_N) * s_w_molar + s_c * s
    MW = x_N * MW_LiBr + (1 - x_N) * MW_H2O
    result = s_molar / MW
    return result
    
def massSpecificHeat(T,x):
    """Inputs:  T = Temperature / [Kelvin]
         x = mass fraction LiBr
Outputs: cp = mass specific heat / [J/kg-K]

Based on Table 6 and equation (3) of reference.
"""
    a = [-14.2094,40.4943,111.135,229.980,1345.26,-0.0141010,0.0124977,-0.000683209]
    m = [2,3,3,3,3,2,1,1]
    n = [0,0,1,2,3,0,3,2]
    t = [0,0,0,0,0,2,3,4]
    Cp_t = 76.0226e3 # [J/kmol-K]
    T_c=647.096 # [K]
    T_0 = 221 # [K] "is a nonlinear parameter of the equations"
    TK = T
    x_N = molefraction(x)
    s=0
    for i in range(len(a)):
        s = s + a[i] * x_N ** m[i] * abs(0.4 - x_N) ** n[i] * (T_c / (TK - T_0)) ** t[i]

    Qu_water = 0.0
    Cp_w_mass = PropsSI('C','T',TK,'Q',Qu_water,'water') # J/kg-K
    Cp_w_molar = Cp_w_mass * MW_H2O
    Cp_molar = (1 - x_N) * Cp_w_molar + Cp_t * s
    MW = x_N * MW_LiBr + (1 - x_N) * MW_H2O
    result = Cp_molar / MW
    return result
    
def twoPhaseProps(h,P,z):
    """Some notes.
    This function returns the quality, temperature and liquid composition of a
    2-phase mixture of liquid lithium bromide-water and water vapor at specific
    enthalpy h, pressure P, and overall composition, z.
    
    h is expected to be in the unit system EES has been configured to operate
    in which can be J/kg or kJ/kg

    P is pressure/bar.

    The compositions, z and x are assumed to be the lithium bromide  mass
    fraction if EES is configured to operate with specific properties on a mass
    basis.  If EES is configured to operate on a  molar basis, z and x are
    assumed to be the lithium bromide mole fraction.
    
    The temperature, T, is assumed to be in the temperature units specified in
    the EES Unit System dialog, which can be C, K, F, or K.

    Q is the quality (or vapor fraction) on a mass basis.
    x is the lithium bromide mass fraction of the liquid phase.

    We observe that all the lithium bromide mass is in the liquid phase.
    Therefore, a mass balance requires (1 - Q) x = z.
    An enthalpy balance gives simply h = (1-Q) h_liquid + Q h_vapor.
    Rearranging, (h - h_liquid) = Q (h_vapor - h_liquid).
    Therefore we have two equations to solve for Q and x.
    Equilibrium requires T_liquid = T_vapor, so we can use exisiting functions.
    """
    P_pascal = P * 1e5
    Q, T, x = 0, 0, 0
    Q = -100	# subcooled
    x = z
    T = temperature(P, x) 
    hL = massSpecificEnthalpy(T,x)
    if (h == hL): Q = 0
    if (h <= hL): return Q, T, x

    Q = 0.1
    for iter in range(100):
        Qlast = Q
        x = z / (1. - Q)
        T = temperature(P, x) 
        hL = massSpecificEnthalpy(T,x)
        hv = 0
        if (h > hL):
            Q_vapor = 1.
            hv = PropsSI('H','T',T,'Q',Q_vapor,'Water')
            hfg = hv - hL
            Q = (h - hL) / (hfg)
            # qq = (x - z) / x
        else:
            Q = 0.
        print("{},h={},P={},z={},Q={},x={},T={},hL={},hv={}"
        .format(iter,h,P,z,Q,x,T,hL,hv))
        if (abs(Q - Qlast) < 0.00001) and (iter > 5):
            break
    return Q, T, x

if __name__ == "__main__":
    # Unit testing
    # Pressure (cf table 6.1, or rather LiBrSS7B.EES (slight difference))
    #for TC,x in ((32.7,0.567),(63.3,0.567),(89.9,0.625),(53.3,0.625),\
    #    (44.7,62.5),(77.0,0.0),(40.2,0.0),(1.3,0.0)):
    #    P = P_LiBrH2O(C2K(TC),x)
    #    print("(T,x) = ({},{}) -> P = {}".format(TC,x,P))
    # Confer documentation for massfraction_LiBrH2O
    for x,w_expected in ((0.1718, 0.50),):
        w = massfraction(x)
        print("Mole fraction {} -> mass fraction {}, expected {}".format(
            x,w,w_expected))
    # Confer documentation for molefraction_LiBrH2O
    for w,x_expected in ((0.5, 0.1718),):
        x = molefraction(w)
        print("Mass fraction {} -> mole fraction {}, expected {}".format(
            w,x,x_expected))
    # Confer documentation for P_LiBrH2O
    for TC,x in ((50.0,0.5),):
        T = C2K(TC)
        P = pressure(T,x)
        P_expected = 0.03486
        print("(T,x) = ({} C, {}) -> P = {} bar".format(TC,x,P))
        print("Expected {}".format(P_expected))
    # Confer documentation for T_LiBrH2O
    for P,x in ((3.5e3 / 1.0e5, 0.5),):
        T = temperature(P,x)
        T_expected = 50.08
        print("(P,x) = {} bar, {} -> T = {} K".format(P,x,K2C(T)))
        print("Expected {} C".format(T_expected))
    # Confer documentation for x_LiBrH2O
    for TC,P in ((50.0, 3.5e3 / 1.0e5),):
        T = C2K(TC)
        x = massFraction(T,P)
        x_expected = 0.4995
        print("(TC,P) = {}, {} -> x = {}".format(TC,P,x))
        print("Expecting {}".format(x_expected))
    # Confer documentation for h_LiBrH2O
    for TC,x,h_expected in ((50,0.5,105e3),):
        T = C2K(TC)
        h = massSpecificEnthalpy(T,x)
        print("T,x = {} C, {} -> h = {} J/kg, expected {}".format(TC,x,h,h_expected))
    for TCvar in np.linspace(20,80,10):
        T = C2K(TCvar)
        h = massSpecificEnthalpy(T,x)
        print("T,x = {:.3f} C, {} -> h = {} J/kg".format(TCvar,x,h))
    TC = 70.
    for xvar in np.linspace(0.2,0.8,10):
        T = C2K(TC)
        h = massSpecificEnthalpy(T,xvar)
        print("T,x = {} C, {:0.3f} -> h = {} J/kg".format(TC,xvar,h))
        
    print("By the way, pure water @ 0 deg C, quality 0 has enthalpy {} J/kg"
        .format(PropsSI('H','T',273.16,'Q',0,'water')))
    T_crit = PropsSI('water','Tcrit')
    P_crit = PropsSI('water','pcrit')
    print("By the way, pure water @ critical point is has {} K, {} Pa"
        .format(T_crit, P_crit))
    h_crit = PropsSI('H','T',T_crit,'P',P_crit,'water')
    print("By the way, pure water @ critical point has enthalpy {} J/kg"
        .format(h_crit))
    
    # Confer documentation for s_LiBrH2O
    for TC,x,s_expected in ((50.0,0.5,0.3519e3),):
        T = C2K(TC)
        s = massSpecificEntropy(T,x)
        print("T,x = {} C, {} -> s = {} J/kg-K, expected {}".format(TC,x,s,s_expected))
    
    # Confer documentation for Cp_LiBrH2O
    for TC,x,Cp_expected in ((50.0,0.5,2.183e3),):
        T = C2K(TC)
        Cp = massSpecificHeat(T,x)
        print("T,x = {} C, {} -> Cp = {} J/kg-K, expected {}"
            .format(TC,x,Cp,Cp_expected))
    
    # Confer documentation for 
    P1, P2 = 10e3/1e5, 1e3/1e5 # [bar]
    TC1 = 70.0 # [C]
    x1 = 0.6 # [kg/kg]
    T1 = C2K(TC1)
    h1 = massSpecificEnthalpy(T1,x1)
    print("T={}, x={} -> h = {}, P (unused) = {}".format(T1,x1,h1,P1))
    # Throttle to lower pressure"
    h2 = h1
    z2 = x1
    QTx2 = twoPhaseProps(h2,P2,z2)
    # Solution:
    QTx2_expected = 0.0146,C2K(48.6),0.6089 # [dim], [C], [dim]
    print("Got QTx = {}, expected {}".format(QTx2, QTx2_expected))
    