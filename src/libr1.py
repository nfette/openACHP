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
T_evap = 1.5
T_cond = 39.9
P_evap = CP.PropsSI('P','T',C2K(T_evap),'Q',1,water)
P_cond = CP.PropsSI('P','T',C2K(T_cond),'Q',1,water)

x2 = 0.624
dx = 0.057
x1 = x2 - dx

T_gen_inlet = libr_props2.Tsat(x1, P_cond, 80)
T_gen_outlet = libr_props2.Tsat(x2, P_cond, T_gen_inlet)
T_abs_inlet_max = libr_props2.Tsat(x2, P_evap, T_gen_outlet)
T_abs_outlet_max = libr_props2.Tsat(x1, P_evap, T_abs_inlet_max)

h_gen_inlet = libr_props2.Hsat(x1, T_gen_inlet)
h_gen_outlet = libr_props2.Hsat(x2,T_gen_outlet)
h_abs_inlet = libr_props2.Hsat(x2,T_abs_inlet_max)
h_abs_outlet = libr_props2.Hsat(x1,T_abs_outlet_max)

# Specify mass flow
m_pump = 0.05
# Mass balance on LiBr
# m_pump * x1 = m_concentrate * x2
m_concentrate = m_pump * x1 / x2
# Mass balance on Water
m_refrig = m_pump - m_concentrate

# Sepcify SHX effectiveness, temperature-wise
Eff_SHX = 0.4
# Compute SHX outlets, assuming concentrate limits heat flow (C_min)
# Neglect pump work for the present.
DeltaT_max = T_gen_outlet - T_abs_outlet_max
DeltaT_SHX_concentrate = Eff_SHX * DeltaT_max
T_SHX_concentrate_outlet = T_gen_outlet - DeltaT_SHX_concentrate
h_SHX_concentrate_outlet \
    = CP.PropsSI('H','T',C2K(T_SHX_concentrate_outlet),'P',P_cond,librname(x2))
Q_SHX = m_concentrate * (h_gen_outlet - h_SHX_concentrate_outlet)
# Expansion valve
h_abs_pre = h_SHX_concentrate_outlet
if h_abs_pre > h_abs_inlet:
    # Pre-cooling is required to reach saturation temperature
    Q_abs_pre_cool = m_concentrate * (h_abs_pre - h_abs_inlet)
    T_abs_pre = np.nan
    # Minimum vapor pressure for absorption to occur
    P_abs_pre = np.inf
else:
    Q_abs_pre_cool = 0
    T_abs_pre = K2C(CP.PropsSI('T','H',h_abs_pre,'P',P_evap,librname(x2)))
    # Minimum vapor pressure for absorption to occur
    P_abs_pre = CP.PropsSI('P','T',C2K(T_abs_pre),'Q',0,librname(x2))

# Heat rejection in absorber: energy balance
h_abs_vapor_inlet = CP.PropsSI('H','P',P_evap,'Q',1,water) - h_w_ref
Q_abs_main = m_refrig * h_abs_vapor_inlet + m_concentrate * h_abs_inlet \
    - m_pump * h_abs_outlet
Q_abs_total = Q_abs_main + Q_abs_pre_cool

# Energy balance in SHX, pump side
D_in = CP.PropsSI('D','T',C2K(T_abs_outlet_max),'Q',0,librname(x1))
DeltaH_pump = (P_cond - P_evap) / D_in
W_pump = m_pump * DeltaH_pump
h_pump_outlet = h_abs_outlet + DeltaH_pump
DeltaH_SHX_pumpside = Q_SHX / m_pump
h_gen_pre = h_pump_outlet + DeltaH_SHX_pumpside
T_gen_pre = K2C(CP.PropsSI('T','P',P_cond,'H',h_gen_pre,librname(x1)))

Q_gen_pre_heat = m_pump * (h_gen_inlet - h_gen_pre)

# Heat input to generator: energy balance
# TODO: this is not saturated maybe?
h_gen_vapor_outlet = CP.PropsSI('H','P',P_cond,'Q',1,water) - h_w_ref
Q_gen_main = m_refrig * h_gen_vapor_outlet + m_concentrate * h_gen_outlet \
    - m_pump * h_gen_inlet
Q_gen_total = Q_gen_main + Q_gen_pre_heat

# Condenser
h_condenser_outlet = CP.PropsSI('H','P',P_cond,'Q',0,water) - h_w_ref
Q_condenser_reject = m_refrig * (h_gen_vapor_outlet - h_condenser_outlet)

# Expansion valve
h_evap_inlet = h_condenser_outlet

# Evaporator
h_evap_outlet = h_abs_vapor_inlet
Q_evap_heat = m_refrig * (h_evap_outlet - h_evap_inlet)

COP = Q_evap_heat / Q_gen_total
ZeroCheck = W_pump + Q_evap_heat + Q_gen_total - Q_condenser_reject - Q_abs_total

if False:
    # Minimum required cooling of concentrate stream
    DeltaH_concentrate = h_abs_inlet - h_gen_outlet
    Q_to_cool_concentrate = m_concentrate * DeltaH_concentrate
    DeltaH_pumpside = h_gen_inlet - h_abs_outlet
    Q_to_raise_dilute = m_pump * DeltaH_pumpside
    # Q_to_raise_dilute is usually larger. So the "actual" generator inlet
    # is subcooled. 
    
    # Temperature-wise effectiveness
    T_SHX_outlet_concentrate = K2C(CP.PropsSI('T','P',P_cond, 'H', h_abs_inlet, librname(x2)))
    
    # The larger of these is the "effectiveness"
    Eff_pumpside = DeltaT_SHX_pumpside / DeltaT_max
    Eff_concentrate = DeltaT_SHX_concentrate / DeltaT_max
    
    DeltaH_gen_preheat = h_gen_inlet - h_gen_pre
    Q_gen_preheat = DeltaH_gen_preheat * m_pump


names = """T_evap T_cond P_evap P_cond
x1 x2
T_gen_inlet T_gen_outlet T_abs_inlet_max T_abs_outlet_max
h_gen_inlet h_gen_outlet h_abs_inlet h_abs_outlet
m_pump m_concentrate m_refrig
Eff_SHX
T_SHX_concentrate_outlet Q_SHX
T_abs_pre h_abs_pre Q_abs_pre_cool P_abs_pre
Q_abs_main Q_abs_total
T_gen_pre
Q_gen_pre_heat Q_gen_main Q_gen_total
Q_condenser_reject Q_evap_heat COP ZeroCheck
string""".split()
vals = [T_evap, T_cond, P_evap, P_cond,
        x1, x2,
        T_gen_inlet, T_gen_outlet, T_abs_inlet_max, T_abs_outlet_max,
        h_gen_inlet, h_gen_outlet, h_abs_inlet, h_abs_outlet,
        m_pump, m_concentrate, m_refrig,
        Eff_SHX,
        T_SHX_concentrate_outlet, Q_SHX,
        T_abs_pre, h_abs_pre, Q_abs_pre_cool, P_abs_pre,
        Q_abs_main, Q_abs_total,
        T_gen_pre,
        Q_gen_pre_heat, Q_gen_main, Q_gen_total,
        Q_condenser_reject, Q_evap_heat, COP, ZeroCheck,
        False]
units = """C C Pa Pa
kg/kg kg/kg
C C C C
J/kg J/kg J/kg J/kg
kg/s kg/s kg/s
K/K
C W
C J/kg W Pa
W W
C
W W W
W W W/W W
none""".split()
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