# -*- coding: utf-8 -*-
"""
Created on Thu Mar 05 18:32:52 2015

@author: nfette
"""
from hw3_8 import *

if __name__ == "__main__":
    T,Qu = C2K(50), 1.0
    g_water = PropsSI("G","T",T,"Q",Qu,"water")
    
    P = 4e3 / 1e5 # [bar]
    x = libr_props.massFraction(T,P)
    g = libr_props.massSpecificGibbs(T,x)
    mu_libr = (1. / x) * (g - (1. - x) * g_water)
    print("""T,P = {} K, {} bar
->
g_water_vapor = {} J/kg,
x = {}, g_liquid = {} J/kg, mu_libr = {} J/kg"""
    .format(T, P, g_water, x,g,mu_libr))
