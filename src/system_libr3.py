# -*- coding: utf-8 -*-
"""
Created on Sat Apr 30 03:26:36 2016

This file solves a thermal system with a chiller.

@author: nfette
"""

import libr3
import HRHX_integral_model

class system1:
    def __init__(self,UA):
        self.UA = UA
        self.chiller = libr3.ChillerLiBr1()
        self.heatStream = HRHX_integral_model.streamExample1(150,-0.1,4179)
        self.absorberRejectStream = HRHX_integral_model.streamExample1(15,1,4179)
        self.condRejectStream = HRHX_integral_model.streamExample1(35,1,4179)
        self.coldStream = HRHX_integral_model.streamExample1(12,-1,4179)
    def calcUA(self):
        self.chiller.iterate1()
        print(self.chiller)
        self.genstream = self.chiller.getGeneratorStream()
        self.genHX = HRHX_integral_model.counterflow_integrator(
            self.genstream,
            self.heatStream)
        self.absStream = self.chiller.getAbsorberStream()
        self.absHX = HRHX_integral_model.counterflow_integrator(
            self.absorberRejectStream,
            self.absStream)
        self.genUA=self.genHX.calcUA(self.chiller.Q_gen_total)
        self.absUA=self.genHX.calcUA(self.chiller.Q_abs_total)
    def optimize(self):
        pass
    #def __repr__(self):
    #    result = self.
        
if __name__ == "__main__":
    sys = system1(100)
    sys.calcUA()
    HRHX_integral_model.plotFlow(sys.genHX)
    HRHX_integral_model.plotFlow(sys.absHX)
    sys.optimize()
    print sys.chiller