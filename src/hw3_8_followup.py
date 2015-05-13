# -*- coding: utf-8 -*-
"""
Created on Wed May 13 10:28:56 2015

@author: nfette
"""

from __future__ import print_function
import ammonia_props
from hw2_1 import CelsiusToKelvin as C2K
import numpy as np
import matplotlib.pyplot as plt
import inspect, os
me = inspect.getfile(inspect.currentframe())
me2 = os.path.basename(me)


if __name__ == "__main__":
    P = 13.5 # bar
    # Average state points on hot and cold sides of SHX from absorp example.
    x1,x2 = 0.3074, 0.38 # ammonia mass fraction
    T1,T2 = 358.1, 326.5
    
    myprops = ammonia_props.AmmoniaProps()
    f1 = myprops.props('TPx')
    
    print("Hot side of SHX:")
    state1 = f1.call(T1, P, x1)
    print(state1)
    dhdT1 = f1.massSpecificHeat(T=T1,P=P,x=x1)
    print("dh/dT = {} kJ/kg-K".format(dhdT1))
    
    print("Cold side of SHX:")
    state2 = f1.call(T2, P, x2)
    print(state2)
    dhdT2 = f1.massSpecificHeat(T=T2,P=P,x=x2)
    print("dh/dT = {} kJ/kg-K".format(dhdT2))
    
    # Parametric plot
    xx = [0.3074, 0.38]
    TT = np.linspace(300,400)
    for x in xx:
        CC = np.empty(len(TT))
        for i in range(len(TT)):
            T = TT[i]
            dhdT = f1.massSpecificHeat(T=T,P=P,x=x)
            CC[i] = dhdT
        print(CC)
        plt.plot(TT, CC,label="x={}".format(x))
    plt.plot([T1],[dhdT1],'o',label="SHX1")
    plt.plot([T2],[dhdT2],'o',label="SHX2")
    plt.ylim([4,5])
    plt.legend(loc='best')
    plt.xlabel('Temperature T / [K]')
    plt.ylabel('Specific heat dh/dT / [kJ/kg-K]')
    plt.savefig('../img/{}.figure{}.png'.format(me2,plt.gcf().number))