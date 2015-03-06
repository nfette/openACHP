# -*- coding: utf-8 -*-
"""
Created on Thu Mar 05 19:14:08 2015

@author: nfette
"""

import ammonia_props
import libr_props
from hw2_1 import CelsiusToKelvin as C2K

def residualsAmmonia(T,P,x,deltaT=1e-4,deltaP=1e-4):
    myprops = ammonia_props.AmmoniaProps()
    f1 = myprops.props(123)
    state1 = f1.call(T,P,x)
    state2 = f1.call(T,P+deltaP,x)
    state3 = f1.call(T+deltaT,P,x)
    dhdp_T = (state2.h - state1.h) / deltaP
    dvdT_p = (state3.v - state1.v) / deltaT
    a = dhdp_T - (state1.v - T * dvdT_p)
    cp = f1.massSpecificHeat(T=T, P=P, x=x)
    dhdT_p = (state3.h - state1.h) / deltaT
    b = dhdT_p - cp
    print(state1)
    print(state2)
    print(state3)
    print("a = {}. b = {}".format(a,b))
    return a,b

def residualsLiBr(T,P,x,deltaT=1e-4,deltaP=1e-4):
    
    state1 = f1.call(T,P,x)
    state2 = f1.call(T,P+deltaP,x)
    state3 = f1.call(T+deltaT,P,x)
    dhdp_T = (state2.h - state1.h) / deltaP
    dvdT_p = (state3.v - state1.v) / deltaT
    a = dhdp_T - (state1.v - T * dvdT_p)
    cp = f1.massSpecificHeat(T=T, P=P, x=x)
    dhdT_p = (state3.h - state1.h) / deltaT
    b = dhdT_p - cp
    print(state1)
    print(state2)
    print(state3)
    print("a = {}. b = {}".format(a,b))
    return a,b


if __name__ == "__main__":
    TC, P, x = 100, 10, 0.5 # C, bar, dim
    # go to 120 C and 11 bar
    residualsAmmonia(C2K(TC), P, x)
    # units are kJ/kg-K
    
    TC, P, x = 100, 1e3 / 1e5, 0.5 # C, bar, dim
    residualsLiBr(C2K(TC), P, x)
    # units are J/kg-K
    