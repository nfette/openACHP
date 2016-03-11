# -*- coding: utf-8 -*-
"""
Created on Fri Mar 04 10:39:51 2016

The purpose of this file is to obtain the saturation properties of aqueous LiBr
from Patel and Klomfar, Int J Refrig, Vol 20, pp 566-578 (2006) directly from
CoolProp.

Note: for enthalpy, a correction is applied.

When inputs are not T,x or P,x, an inverse functions is required, so a
numerical solver is used.

@author: nfette
"""
import CoolProp.CoolProp as CP
from hw2_1 import CelsiusToKelvin as C2K
from hw2_1 import KelvinToCelsius as K2C
from scipy.optimize import fsolve
import numpy as np

librname = lambda(x): 'INCOMP::LiBr[{}]'.format(x)

# TODO: interpolate for guess values

def Tsat(x,P,T_guess=25):
    # CP.PropsSI('T','P',P_evap,'Q',0,libr(x1)) # unsupported inputs
    P_err = lambda(T): CP.PropsSI('P','T',T,'Q',0,librname(x)) - P
    f = fsolve(P_err, C2K(T_guess))
    return K2C(f[0])

TT = np.vectorize(Tsat)
    
def Tsat2(x,H,T_guess=25):
    # CP.PropsSI('T','P',P_evap,'Q',0,libr(x1)) # unsupported inputs
    H_err = lambda(T): CP.PropsSI('H','T',T,'Q',0,librname(x)) - H
    f = fsolve(H_err, C2K(T_guess))
    return f[0]

xglobal=0
def Xsat(T,P,x_guess = 0.4):
    def P_err (x):
        return CP.PropsSI('P','T',C2K(T),'Q',0,librname(x)) - P
    P_err_vec = np.vectorize(P_err)
    f = fsolve(P_err_vec, x_guess)
    return f[0]

Hcorrector = []

def Hsat(x,T):
    """Applies correction based on concentration, then returns saturation
    enthalpy from CoolProp's LiBr-H2O solution. You should check if the
    state is crystalline; if so, the results are not accurate.
    
    Args
    ----
        x (float)
            Concentration of LiBr in liquid, from 0 to 0.75 [kg/kg].
        T (float)
            Equilibrium temperature, from 0 to 227 [C].
    """
    return CP.PropsSI('H','T',C2K(T),'Q',0,librname(x))

def TwoPhase(H, P, x):
    """Properties of a two-phase mixture.
    
    Args
    ----
        H : float
            Specific enthalpy
        P : float
            pressure
        X : float
            Solute concentration per total mass
        
    Returns
    -------
        T : float
            Temperature.
        Q : float
            Vapor quality, mass basis
        Z : float
            Solute concentration in liquid phase"""
    pass

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    P = 101325
    Tees = np.arange(0,30,10)
    exes = np.linspace(0,0.75)
    for T in Tees:
        h = np.nan * exes
        try:
            for i in range(len(exes)):
                x = exes[i]
                h[i] = CP.PropsSI('H','T',C2K(T),'P',P,librname(x))
        except:
            pass
        plt.plot(exes,h,label="{}".format(T))
        
        