# -*- coding: utf-8 -*-
"""
Created on Sun Feb 08 20:46:41 2015

@author: nfette
"""

from CoolProp.CoolProp import PropsSI
import numpy as np

CelsiusConstants = (273.15, 1.0)

def CelsiusToKelvin(TC):
    return TC + CelsiusConstants[0]
    
def KelvinToCelsius(TK):
    return TK - CelsiusConstants[0]
    
def EntropicTemperatureSensibleConstant2(cp, T0, T1):
    Q = cp * (T1 - T0)
    S = cp * (np.log(T1) - np.log(T0))
    T = Q / S
    return T, Q, S

def EntropicTemperatureCondensationAtConstantTP(h_fg, x0, x1, T):
    Q = h_fg * (x1 - x0)
    S = Q / T
    return T, Q, S
    
if __name__ == "__main__":
    fluid = 'water'
    P = 100000 # Pa
    T1 = PropsSI('T', 'P', P, 'Q', 0.5, fluid)
    T0, T2 = 150, 50 # °C
    T0, T2 = map(CelsiusToKelvin, (T0, T2)) # Kelvin
    
    # Approximate (specific heats) method
    T_gas = 0.5 * (T0 + T1)
    Cp_gas = PropsSI('C', 'P', P, 'T', T_gas, fluid)
    T_liquid = 0.5 * (T1 + T2)
    Cp_liquid = PropsSI('C', 'P', P, 'T', T_liquid, fluid)

    x0, x1 = 1, 0
    h_hot = PropsSI('H', 'P', P, 'T', T0, fluid)
    h_f = PropsSI('H', 'P', P, 'Q', 0, fluid)
    h_g = PropsSI('H', 'P', P, 'Q', 1, fluid)
    h_cold = PropsSI('H', 'P', P, 'T', T2, fluid)
    h_desuper = h_hot - h_g
    h_fg = h_g - h_f
    h_subcool = h_f - h_cold
    
    T_de_superheat, Q_de_superheat, S_de_superheat = \
        EntropicTemperatureSensibleConstant2(Cp_gas, T0, T1)
    T_condense, Q_condense, S_condense = \
        EntropicTemperatureCondensationAtConstantTP(h_fg, x0, x1, T1)
    T_subcool, Q_subcool, S_subcool = \
        EntropicTemperatureSensibleConstant2(Cp_liquid, T1, T2)
    
    Q_total = Q_de_superheat + Q_condense + Q_subcool
    S_total = S_de_superheat + S_condense + S_subcool
    T_SA = Q_total / S_total    
    T_SA_Celsius = KelvinToCelsius(T_SA)
    
    print T0, T1, T2
    print h_desuper, h_fg, h_subcool
    print T_gas, Cp_gas    
    print T_liquid, Cp_liquid

    print T_de_superheat, Q_de_superheat, S_de_superheat
    print T_condense, Q_condense, S_condense
    print T_subcool, Q_subcool, S_subcool
    print T_SA, Q_total, S_total
    
    print u"That is, entropic temperature = {} °C".format(T_SA_Celsius)
    