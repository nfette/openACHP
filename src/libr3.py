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
from collections import namedtuple

water = 'HEOS::Water'
librname = lambda(x): 'INCOMP::LiBr[{}]'.format(x)
amm = lambda(x): 'REFPROP::water[{}]&ammonia[{}]'.format(1-x,x)

ProcessPoint = namedtuple("ProcessPoint","fluid m x T P Q H D C")
def nullPP(fluid):
    return ProcessPoint(fluid,1,2,3,None,None,None,None,None)
    
# Units in this file:
# temperature [C]
# enthalpy [J/kg]
# pressure [Pa]
# mass fraction [kg/kg total]
# effectiveness [K/K]

# The LiBr enthalpy is zero at 293.15 K, 1 atm, per
# http://www.coolprop.org/fluid_properties/Incompressibles.html#general-introduction
# We need to evaluate water enthalpy relative to that, but it is not built-in.
# http://www.coolprop.org/coolprop/HighLevelAPI.html#reference-states
T_ref = 20
P_ref = 101325
h_w_ref = CP.PropsSI('H','T',C2K(T_ref),'P',P_ref,water)

class ChillerLiBr1:
    def __init__(self,
                 T_evap=1.5, T_cond=39.9,
                 x1=0.567, x2=0.624,
                 Eff_SHX=0.64, m_pump=0.05):
        """Args
        ----
            T_evap : float
                Evaporator saturation temperature (deg C)
            T_cond : float
                Condenser saturation temperature (deg C)
            x1 : float
                Pump side mass fraction of LiBr in stream (low) [kg/kg]
            x2 : float
                Return side mass fraction of LiBr in stream (high/concentrate)
            Eff_SHX : float
                Effectiveness of the solution heat exchanger (K/K), 0 to 1
            m_pump : float
                Mass flow rate through the solution pump (kg/s)
        """
        self.T_evap = T_evap
        self.T_cond = T_cond
        self.x1 = x1
        self.x2 = x2
        self.m_pump = m_pump
        self.Eff_SHX = Eff_SHX
        self.dx = x1 - x2
        
        self.P_evap = CP.PropsSI('P','T',C2K(T_evap),'Q',1,water)
        self.P_cond = CP.PropsSI('P','T',C2K(T_cond),'Q',1,water)
        
        self.stateLabels = """abs outlet
pump outlet
gen inlet
gen sat. liquid
gen outlet
SHX conc. outlet
abs inlet
abs sat. liquid
gen vapor outlet
cond sat. vapor
cond outlet
evap inlet
evap sat. liquid
evap sat. vapor
evap outlet""".split('\n')
        self.states=dict((k, nullPP('LiBrH2O')) for k in self.stateLabels)
        
        
        self.T_gen_inlet = 0
        self.T_gen_outlet = 0
        self.T_abs_inlet_max = 0
        self.T_abs_outlet_max = 0
        self.h_gen_inlet = 0
        self.h_gen_outlet = 0
        self.h_abs_inlet = 0
        self.h_abs_outlet = 0
        self.m_concentrate = 0
        self.m_refrig = 0
        self.T_SHX_concentrate_outlet = 0
        self.Q_SHX = 0
        self.T_abs_pre = np.nan
        self.h_abs_pre = np.nan
        self.Q_abs_pre_cool = 0
        self.P_abs_pre = np.nan
        self.Q_abs_main = 0
        self.Q_abs_total = 0
        self.T_gen_pre = np.nan
        self.Q_gen_pre_heat = 0
        self.Q_gen_main = 0
        self.Q_gen_total = 0
        self.Q_condenser_reject = 0
        self.Q_evap_heat = 0
        self.COP = 0
        self.W_pump = 0
        self.f = np.inf
        
    def ZeroCheck(self):
        return self.W_pump + self.Q_evap_heat + self.Q_gen_total - self.Q_condenser_reject - self.Q_abs_total
        
    def iterate1(self):
        """Update the internal parameters."""
        self.T_gen_inlet = libr_props2.Tsat(self.x1, self.P_cond, 80)
        self.T_gen_outlet = libr_props2.Tsat(self.x2, self.P_cond,
                                             self.T_gen_inlet)
        self.T_abs_inlet_max = libr_props2.Tsat(self.x2, self.P_evap,
                                                self.T_gen_outlet)
        self.T_abs_outlet_max = libr_props2.Tsat(self.x1, self.P_evap,
                                                 self.T_abs_inlet_max)
        
        self.h_gen_inlet = libr_props2.Hsat(self.x1, self.T_gen_inlet)
        self.h_gen_outlet = libr_props2.Hsat(self.x2, self.T_gen_outlet)
        self.h_abs_inlet = libr_props2.Hsat(self.x2, self.T_abs_inlet_max)
        self.h_abs_outlet = libr_props2.Hsat(self.x1, self.T_abs_outlet_max)
        
        # Mass balance on LiBr
        self.m_concentrate = self.m_pump * self.x1 / self.x2
        # Mass balance on Water
        self.m_refrig = self.m_pump - self.m_concentrate
        self.f = self.m_pump / self.m_refrig

        # Compute SHX outlets, assuming concentrate limits heat flow (C_min)
        # Neglect pump work for the present.
        DeltaT_max = self.T_gen_outlet - self.T_abs_outlet_max
        DeltaT_SHX_concentrate = self.Eff_SHX * DeltaT_max
        self.T_SHX_concentrate_outlet = self.T_gen_outlet \
            - DeltaT_SHX_concentrate
        self.h_SHX_concentrate_outlet = CP.PropsSI('H',
            'T', C2K(self.T_SHX_concentrate_outlet),
            'P', self.P_cond,
            librname(self.x2))
        self.Q_SHX = self.m_concentrate \
            * (self.h_gen_outlet - self.h_SHX_concentrate_outlet)
        
        # Expansion valve
        self.h_abs_pre = self.h_SHX_concentrate_outlet
        if self.h_abs_pre > self.h_abs_inlet:
            # Pre-cooling is required to reach saturation temperature
            self.Q_abs_pre_cool = self.m_concentrate \
                * (self.h_abs_pre - self.h_abs_inlet)
            self.T_abs_pre = np.nan
            # Minimum vapor pressure for absorption to occur
            self.P_abs_pre = np.inf
        else:
            self.Q_abs_pre_cool = 0
            self.T_abs_pre = K2C(CP.PropsSI('T',
                'H', self.h_abs_pre,
                'P', self.P_evap,
                librname(self.x2)))
            # Minimum vapor pressure for absorption to occur
            self.P_abs_pre = CP.PropsSI('P',
                'T', C2K(self.T_abs_pre),
                'Q', 0,
                librname(self.x2))
                
        # Heat rejection in absorber: energy balance
        self.h_abs_vapor_inlet = CP.PropsSI('H',
            'P',self.P_evap,
            'Q',1,
            water) - h_w_ref
        self.Q_abs_main = self.m_refrig * self.h_abs_vapor_inlet \
            + self.m_concentrate * self.h_abs_inlet \
            - self.m_pump * self.h_abs_outlet
        self.Q_abs_total = self.Q_abs_main + self.Q_abs_pre_cool
        
        # Energy balance in SHX, pump side
        D_in = CP.PropsSI('D',
            'T',C2K(self.T_abs_outlet_max),
            'Q',0,
            librname(self.x1))
        DeltaH_pump = (self.P_cond - self.P_evap) / D_in
        self.W_pump = self.m_pump * DeltaH_pump
        self.h_pump_outlet = self.h_abs_outlet + DeltaH_pump
        DeltaH_SHX_pumpside = self.Q_SHX / self.m_pump
        self.h_gen_pre = self.h_pump_outlet + DeltaH_SHX_pumpside
        if self.h_gen_pre > self.h_gen_inlet:
            # Flash steam
            self.T_gen_pre = np.nan
        else:
            self.T_gen_pre = K2C(CP.PropsSI('T',
                'P', self.P_cond,
                'H', self.h_gen_pre,
                librname(self.x1)))
        
        self.Q_gen_pre_heat = self.m_pump * (self.h_gen_inlet - self.h_gen_pre)
        
        # Heat input to generator: energy balance
        self.h_gen_vapor_outlet = CP.PropsSI('H',
            'P', self.P_cond,
            'T', C2K(self.T_gen_inlet), water) - h_w_ref
        self.vapor_superheat = self.T_gen_inlet - self.T_cond
        self.Q_gen_main = self.m_refrig * self.h_gen_vapor_outlet \
            + self.m_concentrate * self.h_gen_outlet \
            - self.m_pump * self.h_gen_inlet
        self.Q_gen_total = self.Q_gen_main + self.Q_gen_pre_heat
        
        # Condenser
        self.h_condenser_outlet = CP.PropsSI('H',
            'P',self.P_cond,
            'Q', 0, water) - h_w_ref
        self.Q_condenser_reject = self.m_refrig * (self.h_gen_vapor_outlet
            - self.h_condenser_outlet)
        
        # Expansion valve
        self.h_evap_inlet = self.h_condenser_outlet
        
        # Evaporator
        self.h_evap_outlet = self.h_abs_vapor_inlet
        self.Q_evap_heat = self.m_refrig * (self.h_evap_outlet
            - self.h_evap_inlet)
        
        self.COP = self.Q_evap_heat / self.Q_gen_total


    def iterate2(self,T_gen_outlet,T_abs_outlet):
        """Resolve the concentrations."""
        pass
    
    def __repr__(self):
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
        Q_condenser_reject Q_evap_heat COP
        W_pump
        self.f
        ZeroCheck
        string""".split()
        vals = [self.T_evap, self.T_cond, self.P_evap, self.P_cond,
            self.x1, self.x2,
            self.T_gen_inlet, self.T_gen_outlet, self.T_abs_inlet_max, self.T_abs_outlet_max,
            self.h_gen_inlet, self.h_gen_outlet, self.h_abs_inlet, self.h_abs_outlet,
            self.m_pump, self.m_concentrate, self.m_refrig,
            self.Eff_SHX,
            self.T_SHX_concentrate_outlet, self.Q_SHX,
            self.T_abs_pre, self.h_abs_pre, self.Q_abs_pre_cool, self.P_abs_pre,
            self.Q_abs_main, self.Q_abs_total,
            self.T_gen_pre,
            self.Q_gen_pre_heat, self.Q_gen_main, self.Q_gen_total,
            self.Q_condenser_reject, self.Q_evap_heat, self.COP,
            self.W_pump,
            self.f,
            self.ZeroCheck(),
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
        W W W/W
        W
        kg/kg
        W
        none""".split()
        return tabulate.tabulate(zip(names,vals,units))
        

if __name__ == "__main__":
    if True:
        # Example 6.1 in the book
        P1,P2 = 673, 7445
        T1 = K2C(CP.PropsSI('T','P',P1,'Q',1,water))
        T2 = K2C(CP.PropsSI('T','P',P2,'Q',1,water))
        c = ChillerLiBr1(T_evap=T1,T_cond=T2)
        c.x2=libr_props2.Xsat(89.9,c.P_cond)
        c.x1=libr_props2.Xsat(32.7,c.P_evap)
        print "Initializing..."
        print c
        print "Iterating..."
        c.iterate1()    
        print c
    if False:
        # Figure 6.3 in the book
        Eff_SHX = np.linspace(0,1)
        COP = np.zeros_like(Eff_SHX)
        for i in range(len(Eff_SHX)):
            c = ChillerLiBr1(Eff_SHX=Eff_SHX[i])
            try:
                c.iterate1()
                COP[i] = c.COP
            except:
                pass
        import matplotlib.pyplot as plt
        plt.plot(Eff_SHX, COP)
        plt.show()
    