# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 21:38:51 2015

@author: nfette
"""
from __future__ import print_function
from hw2_1 import CelsiusToKelvin as C2K
from CoolProp.CoolProp import PropsSI
from numpy import linspace
import libr_props
import matplotlib.pyplot as plt
import inspect, os
me = inspect.getfile(inspect.currentframe())
me2 = os.path.basename(me)

if __name__ == "__main__":
    T,Qu = C2K(50), 1.0
    g_water = PropsSI("G","T",T,"Q",Qu,"water")
    
    mu_1 = []
    g_liquid = []
    xvals = linspace(0.01,0.99)
    for x in xvals:
        g = libr_props.massSpecificGibbs(T,x)
        g_liquid.append(g)
        mu = (1. / x) * (g - (1. - x) * g_water)
        mu_1.append( mu )
        print("T,x = {}, {} -> g = {} J/kg, mu = {}".format(T,x,g,mu))
    
    plt.close('all')
    plt.plot(xvals, mu_1, label="mu")
    plt.plot(xvals, g_liquid, label="g")
    plt.legend(loc='best')
    plt.title("g_water_vapor = {:g} J/kg".format(g_water))
    plt.savefig('../img/{}.figure{}.png'.format(me2,plt.gcf().number))
    