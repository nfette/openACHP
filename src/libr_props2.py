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
    
def Tsat(x,P,T_guess=25):
    # CP.PropsSI('T','P',P_evap,'Q',0,libr(x1)) # unsupported inputs
    P_err = lambda(T): CP.PropsSI('P','T',T,'Q',0,librname(x)) - P
    f = fsolve(P_err, C2K(T_guess))
    return K2C(f[0])

TT = np.vectorize(Tsat)

def Hsat(x,T):
    return CP.PropsSI('H','T',C2K(T),'Q',0,librname(x))