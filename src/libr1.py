# -*- coding: utf-8 -*-
"""
Created on Fri Mar 04 10:04:46 2016

@author: nfette

A single effect LiBr absorption chiller model.
"""

import CoolProp.CoolProp as CP
from hw2_1 import CelsiusToKelvin as C2K
from hw2_1 import KelvinToCelsius as K2C
import libr_props2
import tabulate
import numpy as np
import matplotlib.pyplot as plt

water = 'HEOS::Water'
librname = lambda(x): 'INCOMP::LiBr[{}]'.format(x)
amm = lambda(x): 'REFPROP::water[{}]&ammonia[{}]'.format(1-x,x)

# The LiBr enthalpy is zero at 293.15 K, 1 atm, per
# http://www.coolprop.org/fluid_properties/Incompressibles.html#general-introduction
# We need to evaluate water enthalpy relative to that, but it is not built-in.
# http://www.coolprop.org/coolprop/HighLevelAPI.html#reference-states
T_ref = 20
P_ref = 101325
h_w_ref = CP.PropsSI('H','T',C2K(T_ref),'P',P_ref,water)


# Specify evaporator outlet and condenser outlet temperatures
T_evap = 3
T_cond = 30
P_evap = CP.PropsSI('P','T',C2K(T_evap),'Q',1,water)
P_cond = CP.PropsSI('P','T',C2K(T_cond),'Q',1,water)

x2 = 0.6
dx = 0.1
x1 = x2 - dx

T_gen_inlet = libr_props2.Tsat(x1, P_cond, 80)
T_gen_outlet = libr_props2.Tsat(x2, P_cond, T_gen_inlet)
T_abs_inlet_max = libr_props2.Tsat(x2, P_evap, T_gen_outlet)
T_abs_outlet_max = libr_props2.Tsat(x1, P_evap, T_abs_inlet_max)

# Sepcify SHX effectiveness
Eff_Hx = 0.64
m_pump = 1
# Mass balance on LiBr
# m_pump * x1 = m_concentrate * x2
m_concentrate = m_pump * x1 / x2
# Mass balance on Water
m_refrig = m_pump - m_concentrate
h_gen_outlet = libr_props2.Hsat(x2,T_gen_outlet)
h_abs_inlet = libr_props2.Hsat(x2,T_abs_inlet_max)
DeltaH_concentrate = h_abs_inlet - h_gen_outlet
Q_to_cool_concentrate = m_concentrate * DeltaH_concentrate

h_abs_outlet = libr_props2.Hsat(x2,T_abs_outlet_max)
h_gen_inlet = libr_props2.Hsat(x2,T_gen_inlet)
DeltaH_pumpside = h_gen_inlet - h_abs_outlet
Q_to_raise_dilute = m_pump * DeltaH_pumpside
# Q_to_raise_dilute is usually larger. So the "actual" generator inlet
# is subcooled. Neglect pump work for the present.
DeltaH_SHX_pumpside = -Q_to_cool_concentrate / m_pump
h_gen_pre = h_abs_outlet + DeltaH_SHX_pumpside
T_gen_pre = K2C(CP.PropsSI('T','P',P_cond,'H',h_gen_pre,librname(x1)))

# Temperature-wise effectiveness
T_SHX_outlet_concentrate = K2C(CP.PropsSI('T','P',P_cond, 'H', h_abs_inlet, librname(x2)))
DeltaT_max = T_gen_outlet - T_abs_outlet_max
DeltaT_SHX_pumpside = T_gen_pre - T_abs_outlet_max
DeltaT_SHX_concentrate = T_gen_outlet - T_SHX_outlet_concentrate
# The larger of these is the "effectiveness"
Eff_pumpside = DeltaT_SHX_pumpside / DeltaT_max
Eff_concentrate = DeltaT_SHX_concentrate / DeltaT_max


names = """x1 x2
T_gen_inlet T_gen_outlet T_abs_inlet_max T_abs_outlet_max
m_pump m_concentrate m_refrig
Q_to_cool_concentrate Q_to_raise_dilute
T_gen_pre T_SHX_oulet_concentrate
Eff_pumpside Eff_concentrate""".split()
vals = [x1, x2,
        T_gen_inlet, T_gen_outlet, T_abs_inlet_max, T_abs_outlet_max,
        m_pump, m_concentrate, m_refrig,
        Q_to_cool_concentrate, Q_to_raise_dilute,
        T_gen_pre, T_SHX_outlet_concentrate,
        Eff_pumpside, Eff_concentrate]
units = 'kg/kg kg/kg C C C C kg/s kg/s kg/s W W C C K/K K/K'.split()
print tabulate.tabulate(zip(names,vals,units))




if False:
    fh = np.vectorize(lambda(T): CP.PropsSI('H','T',C2K(T),'P',P_cond,libr_props2.librname(x1)))
    T = np.linspace(60, T_gen_inlet-0.1)
    h = fh(T)
    
    def hTx(T,x):
        return CP.PropsSI('H','T',C2K(T),'Q',0,libr_props2.librname(x))
    x = np.linspace(x1, x2)
    T2 = np.apply_along_axis(libr_props2.TT,0,x,P_cond,80)
    
    fh2 = np.vectorize(hTx)
    h2=fh2(T2,x)
    
    plt.plot(np.append(T,T2),np.append(h,h2))
    plt.show()