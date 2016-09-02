# -*- coding: utf-8 -*-
"""
Created on Fri Feb 06 15:22:52 2015

@author: nfette
"""

CelsiusConstants = (273.15, 1.0)

def CelsiusToKelvin(TC):
    return TC + CelsiusConstants[0]
    
def KelvinToCelsius(TK):
    return TK - CelsiusConstants[0]
    
def COP_cooling_reversible(T_ei, T_ci, T_hi):
    """Inputs: absolute temperatures in ascending order.
    """
    return (T_ei / T_hi) * (T_hi - T_ci) / (T_ci - T_ei)
    
def COP_heating_reversible(T_ei, T_ci, T_hi):
    """Inputs: absolute temperatures in ascending order.
    """
    return 1. + COP_cooling_reversible(T_ei, T_ci, T_hi)
    
def COP_cooling_partial_Tei(T_ei, T_ci, T_hi):
    return (1 / T_hi) * (T_hi - T_ci) / (T_ci - T_ei) \
        + (T_ei / T_hi) * (T_hi - T_ci) / (T_ci - T_ei) ** 2
    
def COP_cooling_partial_Tci(T_ei, T_ci, T_hi):
    return (T_ei / T_hi) * (-1.) / (T_ci - T_ei) \
        - (T_ei / T_hi) * (T_hi - T_ci) / (T_ci - T_ei)**2

def COP_cooling_partial_Thi(T_ei, T_ci, T_hi):
    return (-T_ei / T_hi**2) * (T_hi - T_ci) / (T_ci - T_ei) \
        + (T_ei / T_hi) * (1.) / (T_ci - T_ei)

def main(me,save=False):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib    
    
    plt.close('all')
    matplotlib.rcParams.update({'font.size':11})
    
    # Part a
    plt.figure(1)
    
    Tei, Tci, Thi = 5, 30, 90 # °C
    Tei, Tci, Thi = map(CelsiusToKelvin, (Tei, Tci, Thi))
    COP  = COP_cooling_reversible (Tei,Tci,Thi)
    COPe = COP_cooling_partial_Tei(Tei,Tci,Thi)
    COPc = COP_cooling_partial_Tci(Tei,Tci,Thi)
    COPh = COP_cooling_partial_Thi(Tei,Tci,Thi)
    print(Tei,Tci,Thi, COP, COPe, COPc, COPh)
    
    DeltaT = np.linspace(-60,60,51)
    TTe = Tei + DeltaT
    TTc = Tci + DeltaT
    TTh = Thi + DeltaT
    
    COP_vary_E = COP_cooling_reversible(TTe, Tci, Thi)
    COP_vary_C = COP_cooling_reversible(Tei, TTc, Thi)
    COP_vary_H = COP_cooling_reversible(Tei, Tci, TTh)
        
    plt.semilogy(DeltaT, COP_vary_E, '-o', label="Varying Evaporator (lowest)", markevery=10)
    plt.semilogy(DeltaT, COP_vary_C, '-d', label="Varying Condensor (middle)", markevery=10)
    plt.semilogy(DeltaT, COP_vary_H, '-s', label="Varying Heat input (highest)", markevery=10)
    
    plt.xlabel("Temperature adjustment $\Delta T$ [K]")
    plt.ylabel("Reversible cooling COP")
    plt.legend()
    plt.ylim(0.1,100)

    if save:    
        plt.savefig('../img/{}.figure{}.png'.format(me,plt.gcf().number))
    else:
        plt.show()
    
    plt.figure(2)
    
    DeltaT = np.linspace(-59,59,52)
    TTe = Tei + DeltaT
    TTc = Tci + DeltaT
    TTh = Thi + DeltaT
    COP_vary_E = COP_cooling_reversible(TTe, Tci, Thi)
    COP_vary_C = COP_cooling_reversible(Tei, TTc, Thi)
    COP_vary_H = COP_cooling_reversible(Tei, Tci, TTh)    
    dCOPdE = COP_cooling_partial_Tei(TTe, Tci, Thi)
    dCOPdC = COP_cooling_partial_Tci(Tei, TTc, Thi)
    dCOPdH = COP_cooling_partial_Thi(Tei, Tci, TTh)
    
    DeltaT_midpoints = 0.5 * (DeltaT[1:] + DeltaT[:-1])
    dCOPdE_numerical = np.diff(COP_vary_E) / np.diff(DeltaT)
    dCOPdC_numerical = np.diff(COP_vary_C) / np.diff(DeltaT)
    dCOPdH_numerical = np.diff(COP_vary_H) / np.diff(DeltaT)    
    
    plt.plot(DeltaT, dCOPdE, 'b', label="Varying Evaporator (lowest)", markevery=10)
    plt.plot(DeltaT, dCOPdC, 'g', label="Varying Condensor (middle)", markevery=10)
    plt.plot(DeltaT, dCOPdH, 'r', label="Varying Heat input (highest)", markevery=10)
    plt.plot(DeltaT_midpoints, dCOPdE_numerical, 'o', markevery=10)
    plt.plot(DeltaT_midpoints, dCOPdC_numerical, 'd', markevery=10)
    plt.plot(DeltaT_midpoints, dCOPdH_numerical, 's', markevery=10)
    
    plt.gca().axhline(0,color='k',ls='--')
    plt.xlabel("Temperature adjustment $\Delta T$ [K]")
    plt.ylabel("Reversible cooling COP partial derivative [1/K]")
    plt.legend()
    plt.ylim(-0.2,0.2)
    
    if save:
        plt.savefig('../img/{}.figure{}.png'.format(me,plt.gcf().number))
    else:
        plt.show()
    
if __name__ == "__main__":
    import os
    me = os.path.basename(__file__)
    main(me, save=False)
