# -*- coding: utf-8 -*-
"""
Created on Fri Feb 08 16:48:49 2015

@author: nfette
"""

CelsiusConstants = (273.15, 1.0)

def CelsiusToKelvin(TC):
    return TC + CelsiusConstants[0]
    
def KelvinToCelsius(TK):
    return TK - CelsiusConstants[0]
    
def COP_heat_transformer(T0, T1, T2):
    """Inputs: absolute temperatures in ascending order.
    """
    return (T1 - T0)/ (T1) * (T2) / (T2 - T0)
    
def COP_HT_partial_T0(T0, T1, T2):
    return (-1.) / (T1) * (T2) / (T2 - T0) + (T1 - T0) / T1 * T2 / (T2 - T0)**2
    
def COP_HT_partial_T1(T0, T1, T2):
    return T0 / T1**2 * (T2) / (T2 - T0)

def COP_HT_partial_T2(T0, T1, T2):
    return (T1 - T0) / T1 * (-T0) / (T2 - T0)**2

if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib
    import inspect
    me = inspect.getfile(inspect.currentframe())
    
    plt.close('all')
    matplotlib.rcParams.update({'font.size':11})
    
    # Part a
    plt.figure(1)
    
    T0, T1, T2 = 5, 30, 90 # °C
    T0, T1, T2 = map(CelsiusToKelvin, (T0, T1, T2)) # Kelvin
    COP  = COP_heat_transformer (T0,T1,T2)
    COP0 = COP_HT_partial_T0(T0,T1,T2)
    COP1 = COP_HT_partial_T1(T0,T1,T2)
    COP2 = COP_HT_partial_T2(T0,T1,T2)
    print T0,T1,T2, COP, COP0, COP1, COP2
    
    DeltaT = np.linspace(-60,60,51)
    TT0 = T0 + DeltaT
    TT1 = T1 + DeltaT
    TT2 = T2 + DeltaT
    
    COP_vary_0 = COP_heat_transformer(TT0, T1, T2)
    COP_vary_1 = COP_heat_transformer(T0, TT1, T2)
    COP_vary_2 = COP_heat_transformer(T0, T1, TT2)
    
    plt.plot(DeltaT, COP_vary_0, '-o', label="Varying lowest", markevery=10)
    plt.plot(DeltaT, COP_vary_1, '-d', label="Varying middle", markevery=10)
    plt.plot(DeltaT, COP_vary_2, '-s', label="Varying highest", markevery=10)
    
    plt.xlabel("Temperature adjustment $\Delta T$ [K]")
    plt.ylabel("Reversible heat transformer COP")
    plt.legend()
    plt.ylim(0.1,1)
    #plt.show()
    
    plt.savefig('{}.figure{}.png'.format(me,plt.gcf().number))


    plt.figure(2)
    
    DeltaT = np.linspace(-59,59,52)
    TT0 = T0 + DeltaT
    TT1 = T1 + DeltaT
    TT2 = T2 + DeltaT
    COP_vary_0 = COP_heat_transformer(TT0, T1, T2)
    COP_vary_1 = COP_heat_transformer(T0, TT1, T2)
    COP_vary_2 = COP_heat_transformer(T0, T1, TT2)
    dCOPd0 = COP_HT_partial_T0(TT0, T1, T2)
    dCOPd1 = COP_HT_partial_T1(T0, TT1, T2)
    dCOPd2 = COP_HT_partial_T2(T0, T1, TT2)
    
    DeltaT_midpoints = 0.5 * (DeltaT[1:] + DeltaT[:-1])
    dCOPd0_numerical = np.diff(COP_vary_0) / np.diff(DeltaT)
    dCOPd1_numerical = np.diff(COP_vary_1) / np.diff(DeltaT)
    dCOPd2_numerical = np.diff(COP_vary_2) / np.diff(DeltaT)    
    
    plt.plot(DeltaT, dCOPd0, 'b', label="Varying lowest")
    plt.plot(DeltaT, dCOPd1, 'g', label="Varying middle")
    plt.plot(DeltaT, dCOPd2, 'r', label="Varying highest")
    plt.plot(DeltaT_midpoints, dCOPd0_numerical, 'o', markevery=10)
    plt.plot(DeltaT_midpoints, dCOPd1_numerical, 'd', markevery=10)
    plt.plot(DeltaT_midpoints, dCOPd2_numerical, 's', markevery=10)
    
    plt.gca().axhline(0,color='k',ls='--')
    plt.xlabel("Temperature adjustment $\Delta T$ [K]")
    plt.ylabel("Reversible heat transformer COP partial derivative [1/K]")
    plt.legend()
    #plt.ylim(-0.1,0.1)
    #plt.show()
    
    plt.savefig('{}.figure{}.png'.format(me,plt.gcf().number))
    