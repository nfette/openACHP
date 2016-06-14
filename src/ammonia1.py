# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 13:23:47 2016

A single effect ammonia-water cycle.

@author: nfette
"""

import CoolProp.CoolProp as CP
from hw2_1 import CelsiusToKelvin as C2K
from hw2_1 import KelvinToCelsius as K2C
import tabulate
import numpy as np

from ammonia_props import AmmoniaProps, StateType, State

amm = AmmoniaProps()

class stateIterator(object):
    def __init__(self,chiller):
        self.chiller=chiller
        self.i = iter(chiller.points)
    def __iter__(self):
        return self
    def next(self):
        return self.chiller.__getattribute__(self.i.next())

class AmmoniaChiller(object):
    def __init__(self):
        self.points = """rich abs outlet
rich pump outlet
rich shx outlet
rich gen sat liquid
weak gen outlet
weak shx outlet
weak exp outlet
gen vapor outlet
gen reflux inlet
refrig rect outlet
refrig cond outlet
refrig cehx liquid outlet
refrig exp outlet
refrig evap outlet
refrig cehx sat vapor
refrig cehx vapor outlet
rectifier_liquid
gen_vapor_formation
abs_vapor_final""".replace(" ","_").split('\n')
        vars = """Q_abs,W
Q_gen,W
Q_cond,W
Q_evap,W
W_pump,W
m_rich,kg/s
m_weak,kg/s
m_gen_vapor,kg/s
m_gen_reflux,kg/s
m_refrig,kg/s""".split()
        self.vars = [var.split(',')[0] for var in vars]
        self.units = [var.split(',')[1] for var in vars]
        
        self.rich_abs_outlet = amm.props2(T=400,P=10,x=0.5)
        self.rich_pump_outlet = amm.props2(T=400,P=10,x=0.5)
        self.rich_shx_outlet = amm.props2(T=400,P=10,x=0.5)
        self.rich_gen_sat_liquid = amm.props2(T=400,P=10,x=0.5)
        self.weak_gen_outlet = amm.props2(T=400,P=10,x=0.5)
        self.weak_shx_outlet = amm.props2(T=400,P=10,x=0.5)
        self.weak_exp_outlet = amm.props2(T=400,P=10,x=0.5)
        self.gen_vapor_outlet = amm.props2(T=400,P=10,x=0.5)
        self.gen_reflux_inlet = amm.props2(T=400,P=10,x=0.5)
        self.refrig_rect_outlet = amm.props2(T=400,P=10,x=0.5)
        self.refrig_cond_outlet = amm.props2(T=400,P=10,x=0.5)
        self.refrig_cehx_liquid_outlet = amm.props2(T=400,P=10,x=0.5)
        self.refrig_exp_outlet = amm.props2(T=400,P=10,x=0.5)
        self.refrig_evap_outlet = amm.props2(T=400,P=10,x=0.5)
        self.refrig_cehx_sat_vapor = amm.props2(T=400,P=10,x=0.5)
        self.refrig_cehx_vapor_outlet = amm.props2(T=400,P=10,x=0.5)
        self.rectifier_liquid = amm.props2(T=400,P=10,x=0.5)
        self.gen_vapor_formation = amm.props2(T=400,P=10,x=0.5)
        self.abs_vapor_final = amm.props2(T=400,P=10,x=0.5)
        
        self.Q_abs = 0
        self.Q_gen = 0
        self.Q_cond = 0
        self.Q_evap = 0
        self.W_pump = 0
        self.m_rich = 1
        self.m_weak = 0
        self.m_gen_vapor = 0
        self.m_gen_reflux = 0
        self.m_refrig = 0
        
    def stateTable(self):
        ii = range(len(self.points))
        states = np.zeros_like(ii,dtype=StateType)
        for i,s2 in zip(ii, stateIterator(self)):
            states[i] = s2
        return states
        
    def variablesTable(self):
        dt=[('name','S256'),('unit','S256'),('value','f')]
        thevars = np.zeros_like(self.vars,dtype=dt)
        ii = range(len(thevars))
        for i,var,unit in zip(ii,self.vars,self.units):
            thevars[i] = var,unit,self.__getattribute__(var)
        return thevars
    
    def __repr__(self):
        #return tabulate.tabulate([self.states],StateType.names)
        states = tabulate.tabulate([[p] \
                                    + [s[name] for name in StateType.names]
            for p,s in zip(self.points,self.stateTable())],StateType.names)
        thevars = tabulate.tabulate(self.variablesTable(),
            'name unit value'.split())
        return states + '\n\n' + thevars
    
    def update(self, x_refrig=0.999869,
               T_evap=5.3+273.15,
               T_cond=38.7+273.15,
               Qu_evap=0.998,
               eff_CEHX=0.95,
               T_abs_outlet=37+273.15,
               T_gen_outlet=101+273.15,
               m_rich=0.40,
               eta_pump = 0.8,
               eff_SHX=0.85,
               T_rect = 40.5 + 273.15):
        """Solves the cycle.
        
        Refrigerant cycle step 1:
        1. Input condenser and evaporator outlet temperatures, mass fraction.
        2. Compute condenser and evaporator pressures.
        Refrigerant CEHX:
        1. Input inlet states, effectiveness.
        2. Compute outlet states.
        Refrigerant expansion valve:
        1. Input inlet state
        2. Compute outlet state
        Absorber step 1:
        1. Input absorber pressure and outlet temperature.
        2. Compute outlet mass fraction.
        Generator step 1:
        1. Input generator pressure and outlet temperature.
        2. Compute outlet mass fraction.
        Refrigerant cycle step 2:
        1. Compute mass flow rates.
        Pump:
        1. Input inlet state and outlet pressure.
        2. Compute outlet state.
        SHX:
        1. Input inlet states and effectiveness.
        2. Output outlet states.
        Solution expansion valve:
        1. Input inlet state.
        2. Compute outlet state.
        Absorber step 2:
        1. Input inlet temperature.
        2. Determine heat flow and saturation temperature.
        Generator step 2:
        1. Input inlet temperature.
        2. Determine heat flow and saturation temperature.
        """
        self.refrig_evap_outlet, self.refrig_cond_outlet \
            = self.updateRefrig(x_refrig, T_evap, T_cond, Qu_evap)
        self.refrig_cehx_liquid_outlet, \
            self.refrig_cehx_sat_vapor, \
            self.refrig_cehx_vapor_outlet \
            = self.updateCEHX(self.refrig_cond_outlet,
                              self.refrig_evap_outlet,
                              eff_CEHX)
        self.refrig_exp_outlet \
            = self.updateExpander(self.refrig_cehx_liquid_outlet,
                                  self.refrig_evap_outlet.P)
        self.rich_abs_outlet \
            = self.updateAbsorber1(self.refrig_cehx_vapor_outlet.P,
                                   T_abs_outlet)
        self.weak_gen_outlet \
            = self.updateGenerator1(self.refrig_cond_outlet.P,
                                    T_gen_outlet)
        self.m_rich = m_rich
        x_weak, x_rich = self.weak_gen_outlet.x, self.rich_abs_outlet.x
        self.m_refrig, self.m_weak = self.updateFlowRates(
            m_rich, x_rich, x_weak, x_refrig)
        self.rich_pump_outlet = self.updatePump(self.rich_abs_outlet,
                                                self.refrig_cond_outlet.P,
                                                eta_pump)
        self.rich_shx_outlet, self.weak_shx_outlet = \
            self.updateSHX(self.rich_pump_outlet, self.weak_gen_outlet,
                           eff_SHX,self.m_rich,self.m_weak)
        self.weak_exp_outlet = \
            self.updateExpander(self.weak_shx_outlet,self.refrig_evap_outlet.P)
        # Update some saturated states
        self.rich_gen_sat_liquid = amm.props2(x=self.rich_shx_outlet.x,
                                              P=self.rich_shx_outlet.P,
                                              Qu=0)
        self.gen_vapor_outlet = amm.props2(P=self.rich_gen_sat_liquid.P,
                                           T=self.rich_gen_sat_liquid.T,
                                           Qu=1)
        self.gen_reflux_inlet = amm.props2(P=self.rich_gen_sat_liquid.P,
                                           T=self.rich_gen_sat_liquid.T,
                                           Qu=0)
        self.refrig_rect_outlet = amm.props2(P=self.refrig_cond_outlet.P,
                                             x=x_refrig,
                                             T=T_rect)
        x_vapor, x_liquid = self.gen_vapor_outlet.x, self.gen_reflux_inlet.x
        self.m_gen_vapor,self.m_gen_reflux = \
            self.updateRectifier(self.m_refrig, x_refrig, x_vapor, x_liquid)
        self.rectifier_liquid = amm.props2(T=T_rect,
                                           P=self.refrig_rect_outlet.P,
                                           Qu=0)
        self.gen_vapor_formation = amm.props2(T=T_gen_outlet,
                                              P=self.weak_gen_outlet.P,
                                              Qu=1)
        self.abs_vapor_final = amm.props2(T=self.weak_exp_outlet.T,
                                          P=self.weak_exp_outlet.P,
                                          Qu=1)
        #self.Q_abs = self.m_rich * self. self.m_weak * self.weak_exp_outlet
    def updateRefrig(self, x_refrig, T_evap, T_cond, Qu_evap):
        evap_outlet = amm.props2(T=T_evap,
                                 x=x_refrig,
                                 Qu=Qu_evap)
        cond_outlet = amm.props2(T=T_cond,
                                 x=x_refrig,
                                 Qu=0)
        return evap_outlet, cond_outlet
    def updateCEHX(self, liquid_inlet, vapor_inlet, effectiveness):
        liquid_max = amm.props2(T=vapor_inlet.T,
                                x=liquid_inlet.x,
                                P=liquid_inlet.P)
        vapor_max = amm.props2(T=liquid_inlet.T,
                               x=vapor_inlet.x,
                               P=vapor_inlet.P)
        deltaH_liquid_max = liquid_inlet.h - liquid_max.h
        deltaH_vapor_max = vapor_max.h - vapor_inlet.h
        deltaH_max = min(deltaH_liquid_max, deltaH_vapor_max)
        deltaH = effectiveness * deltaH_max
        liquid_outlet = amm.props2(x=liquid_inlet.x,
                                   P=liquid_inlet.P,
                                   h=liquid_inlet.h - deltaH)
        vapor_outlet = amm.props2(x=vapor_inlet.x,
                                   P=vapor_inlet.P,
                                   h=vapor_inlet.h + deltaH)
        sat = amm.props2(x=vapor_inlet.x,
                         P=vapor_inlet.P,
                         Qu=1)
        return liquid_outlet, sat, vapor_outlet
    def updateExpander(self,inlet,P_outlet):
        outlet = amm.props2(x=inlet.x, h=inlet.h, P=P_outlet)
        return outlet
    def updateAbsorber1(self,P_inlet,T_outlet):
        outlet = amm.props2(P=P_inlet,T=T_outlet,Qu=0)
        return outlet
    def updateGenerator1(self,P_inlet,T_outlet):
        outlet = amm.props2(P=P_inlet,T=T_outlet,Qu=0)
        return outlet
    def updateFlowRates(self,m_rich,x_rich,x_weak,x_refrig):
        # m_rich = m_weak + m_refrig
        # m_rich * x_rich = m_refrig * x_refrig + m_weak * x_weak
        m_weak = m_rich * (x_rich - x_refrig) / (x_weak - x_refrig)
        m_refrig = m_rich * (x_rich - x_weak) / (x_refrig - x_weak)
        return m_refrig, m_weak
    def updatePump(self,inlet,P_outlet,eta_pump):
        """Pump:
        1. Input inlet state and outlet pressure.
        2. Compute outlet state."""
        outlet_ideal = amm.props2(x=inlet.x,P=P_outlet,s=inlet.s)
        deltaH_ideal = outlet_ideal.h - inlet.h
        deltaH = deltaH_ideal / eta_pump
        h_out = inlet.h + deltaH
        outlet = amm.props2(x=inlet.x,P=P_outlet,h=h_out)
        return outlet
    def updateSHX(self,cold_inlet,hot_inlet,effectiveness,m_cold,m_hot):
        """SHX:
        1. Input inlet states and effectiveness.
        2. Output outlet states."""
        cold_outlet = amm.props2(T=hot_inlet.T,
                                x=cold_inlet.x,
                                P=cold_inlet.P)
        hot_outlet = amm.props2(T=cold_inlet.T,
                               x=hot_inlet.x,
                               P=hot_inlet.P)
        deltaH_cold_max = cold_outlet.h - cold_inlet.h
        deltaH_hot_max = hot_inlet.h - hot_outlet.h
        Q_max = min(m_cold * deltaH_cold_max, m_hot * deltaH_hot_max)
        Q = effectiveness * Q_max
        cold_outlet = amm.props2(x=cold_inlet.x,
                                 P=cold_inlet.P,
                                 h=cold_inlet.h + Q / m_cold)
        hot_outlet = amm.props2(x=hot_inlet.x,
                                  P=hot_inlet.P,
                                  h=hot_inlet.h - Q / m_hot)
        # Cold/Hot refer to the streams
        return cold_outlet, hot_outlet
    def updateAbsorber2(self):
        """Absorber step 2:
        1. Input inlet temperature.
        2. Determine heat flow and saturation temperature."""
    def updateGenerator2(self):
        """Generator step 2:
        1. Input inlet temperature.
        2. Determine heat flow and saturation temperature."""
    def updateRectifier(self, m_refrig, x_refrig, x_vapor, x_liquid):
        """Return the mass flow rates of vapor and liquid."""
        # m_refrig = m_vapor - m_liquid
        # m_refrig * x_refrig = m_vapor * x_vapor - m_liquid * x_liquid
        m_vapor = m_refrig * (x_refrig - x_liquid) / (x_vapor - x_liquid)
        m_liquid = m_refrig * (x_refrig - x_vapor) / (x_vapor - x_liquid)
        return m_vapor, m_liquid
    def getPaths(self):
        return [(0,1),
            (1,2),
            (2,3),
            (3,7),
            (7,8),
            (8,3),
            (7,9),
            (9,16),
            (16,8),
            (3,4),
            (4,17),
            (17,7),
            (4,5),
            (5,6),
            (6,0),
            (9,10),
            (10,11),
            (11,12),
            (12,13),
            (13,14),
            (14,15),
            (15,0),
            (15,18),
            (18,6)]
        
if __name__ == "__main__":
    a = AmmoniaChiller()
    print a
    a.update()
    print a

        
