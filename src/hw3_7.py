# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 21:38:51 2015

@author: nfette
"""
from __future__ import print_function
from hw2_1 import CelsiusToKelvin as C2K
from numpy import linspace
import libr_props
import matplotlib.pyplot as plt

if __name__ == "__main__":
    T,Qu = C2K(50), 1.0
    g_water = PropsSI("G","T",T,"Q",Qu,"water")    
    
    mu_1 = []
    xvals = linspace(0,0.99)
    for x in xvals:
        g = libr_props.massSpecificGibbs(T,x)
        mu_1.append( (1. / x) * (g - (1. - x) * g_water) )
        print("T,x = {}, {} -> g = {} J/kg".format(T,x,g))
    
    plt.plot(xvals, mu_1)
    plt.show()
    