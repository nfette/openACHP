# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 21:38:51 2015

@author: nfette
"""
from __future__ import print_function
from hw2_1 import CelsiusToKelvin as C2K
from CoolProp.CoolProp import PropsSI
#from CoolProp import AbstractState, QT_INPUTS
from numpy import linspace
import libr_props
import matplotlib.pyplot as plt
import inspect, os
me = inspect.getfile(inspect.currentframe())
me2 = os.path.basename(me)

def main():
    T,Qu = C2K(50), 1.0
    g_water = PropsSI("G","T",T,"Q",Qu,"water") # mass-specific gibbs energy
    # TODO: identify bug to report to CoolProp and fix it.
    # Maybe just need to expose
    #     double CoolProp::AbstractState::gibbsmolar(void)
    # via
    #     src/wrappers/Python/CoolProp/AbstractState.pxd
    # eg
    #     cpdef double gibbsmolar(self) except *
    # Just to check, we could use this different method, but it has the same problem.
    #water = AbstractState("HEOS","Water")
    #water.update(QT_INPUTS, Qu, T)
    #g_water = water.gibbsmolar() / molarmass

    mu_1 = []
    g_liquid = []
    xvals = []
    for x in linspace(0.01,0.99):
        g = libr_props.massSpecificGibbs(T,x)
        #P = 0
        try:
            P = libr_props.pressure(T, x)
        except:
            print("Crystallized!")
            break
        xvals.append(x)
        g_liquid.append(g)
        mu = (1. / x) * (g - (1. - x) * g_water)
        mu_1.append( mu )
        print("T,x = {}, {} -> P = {} bar, g = {} J/kg, mu = {}".format(T,x,P,g,mu))
    
    plt.close('all')
    plt.xlabel("LiBr mass fraction in liquid phase")
    plt.ylabel("Potential function in liquid phase")
    plt.xlim(0,1)
    plt.plot(xvals, mu_1, label="mu_LiBr")
    plt.plot(xvals, g_liquid, label="g_liquid")
    plt.legend(loc='best')
    plt.title("g_water_vapor = {:g} J/kg".format(g_water))
    plt.savefig('../img/{}.figure{}.png'.format(me2,plt.gcf().number))
    
if __name__ == "__main__":
    main()