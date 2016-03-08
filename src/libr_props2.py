# -*- coding: utf-8 -*-
"""
Created on Fri Mar 04 10:39:51 2016

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

def Xsat(T,P,x_guess = 0.4):
    P_err = lambda(x): CP.PropsSI('P','T',C2K(T),'Q',0,librname(x)) - P
    f = fsolve(P_err, x_guess)
    return f[0]

def Hsat(x,T):
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
    