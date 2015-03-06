# -*- coding: utf-8 -*-
"""
Created on Thu Mar 05 17:26:23 2015

@author: nfette
"""

from __future__ import print_function
import ammonia_props
from hw2_1 import CelsiusToKelvin as C2K

if __name__ == "__main__":
    myprops = ammonia_props.AmmoniaProps()
    f1 = myprops.props('TPx')
    T, P, x = C2K(20.), 10., 0.5 # K, bar, dim
    
    state = f1.call(T, P, x)
    print(state)
    dhdT = f1.massSpecificHeat(T=T,P=P,x=x)
    dhdT_expected = 4.604999985
    print("dh/dT = {}, expected {} kJ/kg-K".format(dhdT,dhdT_expected))
    