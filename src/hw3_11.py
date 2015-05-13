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
    TT = [T, T, T + deltaT]
    PP = [P, P + deltaP, P]
    hh = [0,0,0]
    vv = [0,0,0]
    for i in range(3):
        # TODO: make this a function
        h_liquid = libr_props.massSpecificEnthalpy(T,x)
        v_liquid = 1. / libr_props.massDensity(T,x) # m3/kg
        x_liquid = libr_props.massFraction(T,P)
        cp_liquid = libr_props.massSpecificHeat(T,x)
        # Now let's figure out the other thing
        Qu = 1. - x / x_liquid
        Q_vapor = 1
        h_vap = PropsSI('H','T',T,'Q',Q_vapor,'Water')
        v_vap = PropsSI('D','T',T,'Q',Q_vapor,'Water')
        cp_vap = PropsSI('Cp','T',T,'Q',Q_vapor,'Water')
        h = (1.0 - Qu) * h_liquid + Qu * h_vap
        v = (1.0 - Qu) * v_liquid + Qu * v_vap
        hh.append(h)
        vv.append(v)
    dhdp_T = (h[1] - h[0]) / deltaP
    dvdT_p = (v[2] - v[0]) / deltaT
    a = dhdp_T - (v[0] - T * dvdT_p)    
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
    