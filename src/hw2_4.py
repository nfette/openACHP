# -*- coding: utf-8 -*-
"""
Created on Sun Feb 08 17:34:42 2015

@author: nfette
"""

CelsiusConstants = (273.15, 1.0)

def CelsiusToKelvin(TC):
    return TC + CelsiusConstants[0]
    
def KelvinToCelsius(TK):
    return TK - CelsiusConstants[0]
    
def EntropicTemperatureSensibleConstant(T0, T1):
    DeltaT = T1 - T0
    denum = np.log(T1) - np.log(T0)
    return DeltaT / denum

if __name__ == "__main__":
    import numpy as np
    import inspect
    me = inspect.getfile(inspect.currentframe())
    
    T0, T1 = 50, 30 # °C
    T0, T1 = map(CelsiusToKelvin, (T0, T1)) # Kelvin
    TSA = EntropicTemperatureSensibleConstant(T0, T1)
    TSA_celsius = KelvinToCelsius(TSA)
    print T0, T1, TSA, TSA_celsius
    
    