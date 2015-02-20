# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 14:42:55 2015

@author: nfette
"""
from __future__ import print_function
import ammonia_props
from hw2_1 import CelsiusToKelvin as C2K


if __name__ == "__main__":
    print("(a)")
    myprops = ammonia_props.AmmoniaProps()
    f1 = myprops.props('TPx')
    T, P, x = C2K(50.), 10., 0.0 # K, bar, dim    
    state = f1.call(T, P, x)    
    print(state)
    print(f1.getOutputUnits())
    
    print("(b)")
    Meff = state.molarMass()
    hbar = state.h * Meff
    print("Meff = {} kg/kmol".format(Meff))
    print("hbar = {} kJ/kmol".format(hbar))
    
    print("(c)")
    dhdx = f1.dhdxetc(T=T,P=P,x=x)
    print(dhdx)
    print("kJ/kg")
    
    print("(d)")
    x = state.x
    h_ideal1 = f1.call(T,P,1.0).h
    h_ideal2 = f1.call(T,P,0.0).h
    h_ideal = x * h_ideal1 + (1 - x) * h_ideal2
    h_mix = state.h - h_ideal
    
    print("(e)")
    print("h1bar = {}".format(dhdx.h1 * ammonia_props.molecular_mass_ammonia))
    print("h1bar = {}".format(dhdx.h2 * ammonia_props.molecular_mass_water))
    
    print("(f)")
    print("h_mix_bar = {}".format(h_mix * Meff))
    print("kJ/kmol")
    