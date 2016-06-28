# -*- coding: utf-8 -*-
"""
Created on Sat Apr 30 03:26:36 2016

This file solves a thermal system with a chiller.

@author: nfette
"""

import libr3
import HRHX_integral_model
import tabulate
import CoolProp.CoolProp as CP
from hw2_1 import CelsiusToKelvin as C2K
from hw2_1 import KelvinToCelsius as K2C

class system1:
    def __init__(self,UA):
        self.UA = UA
        self.chiller = libr3.ChillerLiBr1()
        self.heatStream = HRHX_integral_model.streamExample1(120,0.3,4179)
        self.absorberRejectStream = HRHX_integral_model.streamExample1(15,1,4179)
        self.condRejectStream = HRHX_integral_model.streamExample1(35,5,4179)
        self.coldStream = HRHX_integral_model.streamExample1(12,1,4179)
    def calcUA(self):
        self.chiller.iterate1()
        print(self.chiller)
                
        self.coldStream = self.makeColdStream(self.chiller.Q_evap_heat)
        self.absorberRejectStream = self.makeRejectStream(32, 5, self.chiller.Q_abs_total)
        self.condRejectStream = self.makeRejectStream(32,5, self.chiller.Q_condenser_reject)
        
        self.genstream = self.chiller.getGeneratorStream()
        self.genHX = HRHX_integral_model.counterflow_integrator(self.genstream,self.heatStream,useHotT=True)
        self.absStream = self.chiller.getAbsorberStream()
        self.absHX = HRHX_integral_model.counterflow_integrator(self.absorberRejectStream,self.absStream)
        self.condStream = self.chiller.getCondenserStream()
        self.condHX = HRHX_integral_model.counterflow_integrator(self.condRejectStream,self.condStream,useHotT=True)
        self.condHX.Qmax = self.chiller.Q_condenser_reject
        self.evapStream = self.chiller.getEvaporatorStream()
        self.evapHX = HRHX_integral_model.counterflow_integrator(self.evapStream,self.coldStream)
        
        self.genUA=self.genHX.calcUA(self.chiller.Q_gen_total)
        self.absUA=self.absHX.calcUA(self.chiller.Q_abs_total)
        self.condUA=self.condHX.calcUA(self.chiller.Q_condenser_reject)
        self.evapUA=self.evapHX.calcUA(self.chiller.Q_evap_heat)
    
    def makeColdStream(self,Q):
        # Flow rate for evaporator is specified by inlet and outlet temperature...
        T_in,T_out = 12, 7
        cp = CP.PropsSI('C','T',C2K(T_in),'Q',0,'water')
        print "cp is ",cp
        dh = cp * (T_out - T_in)
        m = -Q / dh
        return HRHX_integral_model.streamExample1(T_in,m,cp)
        
    def makeRejectStream(self,T_in,dT,Q):
        # Flow rate for reject streams is specified by inlet and outlet temperatures...
        cp = CP.PropsSI('C','T',C2K(T_in),'Q',0,'water')
        print "cp is ",cp
        dh = cp * dT
        m = Q / dh
        return HRHX_integral_model.streamExample1(T_in,m,cp)
    
    def optimize(self):
        pass
    
    def __repr__(self):
        names = """genQ
        genQmax
        genUA
        absQ
        absQmax
        absUA
        condQ
        condQmax
        condUA
        evapQ
        evapQmax
        evapUA""".split()
        vals = [self.chiller.Q_gen_total,
                self.genHX.Qmax,
                self.genUA,
                self.chiller.Q_abs_total,
                self.absHX.Qmax,
                self.absUA,
                self.chiller.Q_condenser_reject,
                self.condHX.Qmax,
                self.condUA,
                self.chiller.Q_evap_heat,
                self.evapHX.Qmax,
                self.evapUA]
        units = "W W W/K".split() * 4
        return tabulate.tabulate(zip(names,vals,units))
        
if __name__ == "__main__":
    sys = system1(100)
    sys.calcUA()
    HRHX_integral_model.plotFlow(sys.genHX, None, sys.chiller.Q_gen_total)
    plt.title("Generator")
    HRHX_integral_model.plotFlow(sys.absHX, None, sys.chiller.Q_abs_total)
    plt.title("Absorber")
    HRHX_integral_model.plotFlow(sys.condHX, None, sys.chiller.Q_condenser_reject)
    plt.title("Condenser")
    HRHX_integral_model.plotFlow(sys.evapHX, None, sys.chiller.Q_evap_heat)
    plt.title("Evaporator")
    sys.optimize()
    print sys