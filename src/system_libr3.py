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
import scipy.optimize

def makeSystem():
        return system1(chiller=libr3.ChillerLiBr1(),
        heatStream=HRHX_integral_model.streamExample1(120,0.3,4179),
        absorberRejectStream=HRHX_integral_model.streamExample1(15,1,4179),
        condRejectStream=HRHX_integral_model.streamExample1(35,1,4179),
        coldStream=HRHX_integral_model.streamExample1(12,1,4179))

class system1:
    def __init__(self,chiller,heatStream,absorberRejectStream,condRejectStream,coldStream):
        self.chiller = chiller
        self.heatStream = heatStream
        self.absorberRejectStream = absorberRejectStream
        self.condRejectStream = condRejectStream
        self.coldStream = coldStream
        
    def calcUA(self,x):
        
        if x is not None:
            mdot,T_evap,T_cond,x1,x2 = x
            self.chiller.m_pump=mdot
            self.chiller.setT_cond(T_cond)
            self.chiller.setT_evap(T_evap)
            self.chiller.x1 = x1
            self.chiller.x2 = x2
        
        self.chiller.iterate1()
        #print(self.chiller)
                
        self.coldStream = self.makeColdStream(self.chiller.Q_evap_heat)
        self.absorberRejectStream = self.makeRejectStream(32, 5, self.chiller.Q_abs_total)
        self.condRejectStream = self.makeRejectStream(20, 5, self.chiller.Q_condenser_reject)        
        
        self.genstream = self.chiller.getGeneratorStream()
        self.genHX = HRHX_integral_model.counterflow_integrator(self.genstream,self.heatStream,useHotT=True)

        self.absStream = self.chiller.getAbsorberStream()        
        self.absHX = HRHX_integral_model.counterflow_integrator(self.absorberRejectStream,self.absStream)
        
        self.condStream = self.chiller.getCondenserStream()        
        self.condHX = HRHX_integral_model.counterflow_integrator(self.condRejectStream,self.condStream,useHotT=True)
        #self.condHX.Qmax = self.chiller.Q_condenser_reject
        
        self.evapStream = self.chiller.getEvaporatorStream()        
        self.evapHX = HRHX_integral_model.counterflow_integrator(self.evapStream,self.coldStream)
        
        print "Generator..."
        self.genUA,self.genEff = self.genHX.calcUA(self.chiller.Q_gen_total,True)
        print "Absorber..."
        self.absUA,self.absEff = self.absHX.calcUA(self.chiller.Q_abs_total,True)
        print "Condenser..."
        self.condUA,self.condEff = self.condHX.calcUA(self.chiller.Q_condenser_reject,True)
        print "Evaporator..."
        self.evapUA,self.evapEff = self.evapHX.calcUA(self.chiller.Q_evap_heat,True)
        
        total = self.genUA + self.absUA + self.condUA + self.evapUA
        return total
    
    def makeColdStream(self,Q):
        # Flow rate for evaporator is specified by inlet and outlet temperature...
        T_in,T_out = 12, 7
        cp = CP.PropsSI('C','T',C2K(T_in),'Q',0,'water')
        #print "cp is ",cp
        dh = cp * (T_out - T_in)
        m = -Q / dh
        return HRHX_integral_model.streamExample1(T_in,m,cp)
        
    def makeRejectStream(self,T_in,dT,Q):
        # Flow rate for reject streams is specified by inlet and outlet temperatures...
        cp = CP.PropsSI('C','T',C2K(T_in),'Q',0,'water')
        #print "cp is ",cp
        dh = cp * dT
        m = Q / dh
        return HRHX_integral_model.streamExample1(T_in,m,cp)
    
    def objective(self,x):
        print "Objective: x = {}".format(x)
        mdot,T_evap,T_cond,x1,x2 = x
        try:
            UA=self.calcUA(x)
            Q = self.chiller.Q_evap_heat
            print " -> Q = {}, UA = {}".format(Q,UA)
            return -Q
        except:
            print "An error occurred in calcUA(). Reporting zeros."
            return 0
                  
    def constraint1(self,x,UAgoal):
        # Result is zero iff: UA is the given value.
        print "Constraint1: x={}, UAgoal={}".format(x,UAgoal)
        try:
            UA = self.calcUA(x)
            print "UA at x = {}".format(UA)
            return UAgoal - UA
        except:
            print "An error occurred in calcUA(). Reporting zeros."
            return 0

    def constraint2(self,x):
        # Result is >=0 iff: x2 > x1
        mdot,T_evap,T_cond,x1,x2 = x
        return x2 - x1
        
    def constraint3(self,x):
        # Result is >=0 iff: T_cond > T_evap
        mdot,T_evap,T_cond,x1,x2 = x
        return T_cond - T_evap

    def constraint4(self,x):
        # Result is >=0 iff: mdot > 0
        mdot,T_evap,T_cond,x1,x2 = x
        return mdot
        
    def constraint5(self,x):
        # Result is >=0 iff: x1 > 0
        mdot,T_evap,T_cond,x1,x2 = x
        return x1
        
    def constraint6(self,x):
        # Result is >=0 iff: x2 < 1
        mdot,T_evap,T_cond,x1,x2 = x
        return 0.7 - x2
    
    def constraint7(self,x):
        # Result is >=0 iff: T_evap > 0
        mdot,T_evap,T_cond,x1,x2 = x
        return T_evap
        
    def constraint8(self,x):
        # Result is >=0 iff: T_evap > 0
        mdot,T_evap,T_cond,x1,x2 = x
        return 65 - T_cond
        
    def constraint9(self,x):
        # Result is >=0 iff: T_evap > 0
        mdot,T_evap,T_cond,x1,x2 = x
        return 0.3 - mdot
        
    def optimize(self):
        # Minimize total UA given input streams
        m_pump,T_evap,T_cond,x1,x2 = 0.05, 1.5, 39.9, 0.567, 0.624
        m_pump,T_evap,T_cond,x1,x2 = 0.06213797,1.49332455,39.90638902,0.61209421,0.6791503
        #Eff_SHX=0.64
        
        x0 = m_pump, T_evap, T_cond, x1, x2        
        UAgoal = 5037.68247933
        
        constraints=[dict(type='ineq',fun=self.constraint1,args=(UAgoal,)),
                     dict(type='ineq',fun=self.constraint2,jac=lambda(x):[0,0,0,-1,1]),
                     dict(type='ineq',fun=self.constraint3,jac=lambda(x):[0,-1,1,0,0]),
                     dict(type='ineq',fun=self.constraint4,jac=lambda(x):[1,0,0,0,0]),
                     dict(type='ineq',fun=self.constraint5,jac=lambda(x):[0,0,0,1,0]),
                     dict(type='ineq',fun=self.constraint6,jac=lambda(x):[0,0,0,0,-1]),
                     dict(type='ineq',fun=self.constraint7,jac=lambda(x):[0,1,0,0,0]),
                     dict(type='ineq',fun=self.constraint8,jac=lambda(x):[0,0,-1,0,0]),
                     dict(type='ineq',fun=self.constraint9,jac=lambda(x):[-1,0,0,0,0])
                    ]
                
        opt = scipy.optimize.minimize(self.objective,x0,constraints=constraints,
                                      method='COBYLA',
                                      options=dict(disp=True,maxiter=10))
        return opt
    
    def __repr__(self):
        names = """genQ
        genEff
        genUA
        absQ
        absEff
        absUA
        condQ
        condEff
        condUA
        evapQ
        evapEff
        evapUA""".split()
        vals = [self.chiller.Q_gen_total,
                self.genEff,
                self.genUA,
                self.chiller.Q_abs_total,
                self.absEff,
                self.absUA,
                self.chiller.Q_condenser_reject,
                self.condEff,
                self.condUA,
                self.chiller.Q_evap_heat,
                self.evapEff,
                self.evapUA]
        units = "W W/W W/K".split() * 4
        return tabulate.tabulate(zip(names,vals,units))
    
    def display(sys):
        import matplotlib.pyplot as plt
        HRHX_integral_model.plotFlow(sys.genHX, None, sys.chiller.Q_gen_total)
        plt.title("Generator")
        HRHX_integral_model.plotFlow(sys.absHX, None, sys.chiller.Q_abs_total)
        plt.title("Absorber")
        HRHX_integral_model.plotFlow(sys.condHX, None, sys.chiller.Q_condenser_reject)
        plt.title("Condenser")
        HRHX_integral_model.plotFlow(sys.evapHX, None, sys.chiller.Q_evap_heat)
        plt.title("Evaporator")    
        print sys
        
def main():

    sys = makeSystem()

    import matplotlib.pyplot as plt
    import numpy as np
    from exceptions import ValueError
    
    if False:
        fig = plt.figure()
        plt.title("Generator")
        for mdot in np.linspace(0.01,0.5):
            try:
                sys.chiller.m_pump = mdot
                print sys.calcUA(None)
                HRHX_integral_model.plotFlow(sys.genHX, fig, sys.chiller.Q_gen_total)
            except ValueError as e:
                print e
                break

    #try:
    m_pump,T_evap,T_cond,x1,x2 = 0.06213797,1.49332455,39.90638902,0.61209421,0.6791503
    x0 = m_pump, T_evap, T_cond, x1, x2
    print sys.calcUA(x0)
    sys.display()
    opt = sys.optimize()
    #finally:
    sys.display()
    
if __name__ == "__main__":
    main()