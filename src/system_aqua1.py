# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 23:39:27 2016

@author: nfette
"""

import ammonia1
import HRHX_integral_model

def makeSystem():
    # Units: K, kg/s, kW/kg-K
    return system_aqua1(chiller=ammonia1.AmmoniaChiller(),
    heatStream=HRHX_integral_model.streamExample1(400,1,4.179),
    absorberRejectStream=HRHX_integral_model.streamExample1(305,3,4.179),
    condRejectStream=HRHX_integral_model.streamExample1(305,5,4.179),
    coldStream=HRHX_integral_model.streamExample1(285,4,4.179),
    rectifierRejectStream=HRHX_integral_model.streamExample1(305,0.15,4.179))

class system_aqua1:
    def __init__(self,chiller,heatStream,absorberRejectStream,condRejectStream,
                 coldStream,rectifierRejectStream):
        self.chiller = chiller
        self.heatStream = heatStream
        self.absorberRejectStream = absorberRejectStream
        self.condRejectStream = condRejectStream
        self.coldStream = coldStream
        self.rectifierRejectStream = rectifierRejectStream

        self.genEff = 0
        self.genUA = 0
        
        self.absEff = 0
        self.absUA = 0
        
        self.condEff = 0
        self.condUA = 0
        
        self.evapEff = 0
        self.evapUA = 0
        
        self.rectEff = 0
        self.rectUA = 0
        
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
        evapUA
        rectQ
        rectEff
        rectUA""".split()
        vals = [self.chiller.Q_gen,
                self.genEff,
                self.genUA,
                self.chiller.Q_abs,
                self.absEff,
                self.absUA,
                self.chiller.Q_cond,
                self.condEff,
                self.condUA,
                self.chiller.Q_evap,
                self.evapEff,
                self.evapUA,
                self.chiller.Q_reflux,
                self.rectEff,
                self.rectUA]
        units = "W W/W W/K".split() * 5
        return tabulate.tabulate(zip(names,vals,units))

    def calcUA(self):
        #self.coldStream = self.makeColdStream(self.chiller.Q_evap_heat)
        #self.absorberRejectStream = self.makeRejectStream(32, 5, self.chiller.Q_abs_total)
        #self.condRejectStream = self.makeRejectStream(20, 5, self.chiller.Q_condenser_reject)        
        
        self.genstream = self.chiller.getGeneratorStream()
        self.genHX = HRHX_integral_model.counterflow_integrator(
            self.genstream,self.heatStream,useHotT=True)

        self.absStream = self.chiller.getAbsorberStream()        
        self.absHX = HRHX_integral_model.counterflow_integrator(
            self.absorberRejectStream,self.absStream)
        
        self.condStream = self.chiller.getCondenserStream()        
        self.condHX = HRHX_integral_model.counterflow_integrator(
            self.condRejectStream,self.condStream,useHotT=True)
        #self.condHX.Qmax = self.chiller.Q_condenser_reject
        
        self.evapStream = self.chiller.getEvaporatorStream()        
        self.evapHX = HRHX_integral_model.counterflow_integrator(
            self.evapStream,self.coldStream)
        
        self.rectifierStream = self.chiller.getRectifierStream()
        self.rectHX = HRHX_integral_model.counterflow_integrator(
            self.rectifierRejectStream,self.rectifierStream)
        
        print "Generator..."
        self.genUA,self.genEff = self.genHX.calcUA(self.chiller.Q_gen,True)
        print "Absorber..."
        self.absUA,self.absEff = self.absHX.calcUA(-self.chiller.Q_abs,True)
        print "Condenser..."
        self.condUA,self.condEff = self.condHX.calcUA(-self.chiller.Q_cond,True)
        print "Evaporator..."
        self.evapUA,self.evapEff = self.evapHX.calcUA(self.chiller.Q_evap,True)
        print "Rectifier..."
        self.rectUA,self.rectEff = self.rectHX.calcUA(-self.chiller.Q_reflux,True)
        
        total = self.genUA + self.absUA + self.condUA + self.evapUA + self.rectUA
        return total
        
    def display(sys):
        HRHX_integral_model.plotFlow(sys.genHX, None, sys.chiller.Q_gen)
        plt.title("Generator")
        HRHX_integral_model.plotFlow(sys.absHX, None, -sys.chiller.Q_abs)
        plt.title("Absorber")
        HRHX_integral_model.plotFlow(sys.condHX, None, -sys.chiller.Q_cond)
        plt.title("Condenser")
        HRHX_integral_model.plotFlow(sys.evapHX, None, sys.chiller.Q_evap)
        plt.title("Evaporator")
        HRHX_integral_model.plotFlow(sys.rectHX, None, -sys.chiller.Q_reflux)
        plt.title("Rectifier")
        
        Q = [sys.chiller.Q_gen,sys.chiller.Q_abs,sys.chiller.Q_cond,
              sys.chiller.Q_evap,sys.chiller.Q_reflux]
        UA = [sys.genUA,sys.absUA,sys.condUA,sys.evapUA,sys.rectUA]
        component = "Generator Absorber Condenser Evaporator Rectifier".split()
        width = 1
        
        plt.figure()
        bar1=plt.bar(range(5),Q,width,tick_label=component)
        plt.figure()
        bar2=plt.bar(range(5),UA,width,tick_label=component)
        print sys
        
if __name__ == "__main__":
    sys = makeSystem()
    sys.chiller.update()
    print sys
    sys.calcUA()
    sys.display()
