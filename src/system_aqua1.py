# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 23:39:27 2016

@author: nfette
"""

import scipy
import tabulate
from HRHX_integral_model import plotFlow as plotFlow, \
    counterflow_integrator as hxci, \
    streamExample1 as se1
import ammonia1

def makeBoundary(x):
    t0,m0,t1,m1,t2,m2,t3,m3,t4,m4 = x
    # Units: K, kg/s, kW/kg-K
    cp = 4.179
    return boundary(
    heat=se1(t0,m0,cp),
    absorberReject=se1(t1,m1,cp),
    condReject=se1(t2,m2,cp),
    cold=se1(t3,m3,cp),
    rectifierReject=se1(t4,m4,cp))

def makeChiller(x):
    chiller = ammonia1.AmmoniaChiller()
    m_rich, T_evap, T_cond, T_rect, T_abs_outlet, T_gen_outlet = x        
    chiller.update(m_rich = m_rich,
                   T_evap = T_evap,
                   T_cond = T_cond,
                   T_rect = T_rect,
                   T_abs_outlet = T_abs_outlet,
                   T_gen_outlet = T_gen_outlet)
    return chiller

class boundary(object):
    def __init__(self,heat,absorberReject,condReject,
                 cold,rectifierReject):
        self.heat            = heat
        self.absorberReject  = absorberReject
        self.condReject      = condReject
        self.cold            = cold
        self.rectifierReject = rectifierReject
    def __repr__(self):
        import tabulate
        result = []
        for name in "heat absorberReject condReject cold rectifierReject".split():
            stream = self.__getattribute__(name)
            result.append((name, stream.T_inlet, stream.mdot, stream.cp))
        return tabulate.tabulate(result,"stream T_inlet mdot cp".split())

class system(object):
    def __init__(self,boundary,chiller):
        self.boundary = boundary
        self.chiller = chiller
        
        self.hxs = {}
        self.hxs['gen'] = hxci(chiller.getGeneratorStream(),boundary.heat)
        self.hxs['rect'] = hxci(boundary.rectifierReject, chiller.getRectifierStream())
        self.hxs['abs'] = hxci(boundary.absorberReject, chiller.getAbsorberStream())
        self.hxs['cond'] = hxci(boundary.condReject, chiller.getCondenserStream())
        self.hxs['evap'] = hxci(chiller.getEvaporatorStream(), boundary.cold)
        self.Q = {'gen':chiller.Q_gen,
                  'rect':-chiller.Q_reflux,
                  'abs':-chiller.Q_abs,
                  'cond':-chiller.Q_cond,
                  'evap':chiller.Q_evap}
        
        self.totalUA = 0
        self.data = []
        for name in self.hxs:
            self.hxs[name].calcQmax()
            deltaT, epsilon, UA = self.hxs[name].calcUA2(self.Q[name])
            self.data.append((name, deltaT, epsilon, UA, self.Q[name]))
            self.totalUA += UA
    
    def __repr__(self):
        result = "{}\ntotalUA = {}".format(
            tabulate.tabulate(self.data, "name deltaT epsilon UA Q".split()),
            self.totalUA)
        return result
        
    def display(self):
        import matplotlib.pyplot as plt
        for name in self.hxs:
            plotFlow(self.hxs[name], None, self.Q[name])
            plt.title(name)
        
#        Q = [sys.chiller.Q_gen,sys.chiller.Q_abs,sys.chiller.Q_cond,
#              sys.chiller.Q_evap,sys.chiller.Q_reflux]
#        UA = [sys.genUA,sys.absUA,sys.condUA,sys.evapUA,sys.rectUA]
#        component = "Generator Absorber Condenser Evaporator Rectifier".split()
#        width = 1
#        
#        plt.figure()
#        bar1=plt.bar(range(5),Q,width,tick_label=component)
#        plt.figure()
#        bar2=plt.bar(range(5),UA,width,tick_label=component)

def main():
    # Boundary
    xB0 = [400,1,
      305,3,
      305,5,
      285,4,
      305,0.15]
    bdry = makeBoundary(xB0)
    print bdry
    # Chiller
    xC0 = [0.40, 278.45, 311.85, 313.65, 310.15, 374.15]
    ch = makeChiller(xC0)
    print ch
    # System
    sys = system(bdry,ch)
    print sys
    sys.display()
    
if __name__ == "__main__":
    main()