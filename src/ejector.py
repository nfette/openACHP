# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 16:41:33 2016

@author: nfette
"""

import CoolProp as CP
import tabulate
import numpy as np
from hw2_1 import CelsiusToKelvin as C2K, KelvinToCelsius as K2C
import HRHX_integral_model
import matplotlib.pyplot as plt

def safe_keyed_output(state,key):
    try:
        return state.keyed_output(key)
    except:
        return np.nan

class EjectorCycle(object):
    def __init__(self,T_hot,T_reject,T_cold,m_dot,
                 eta_nozzle=0.85, eta_mixer=0.70, eta_diffuser=0.70,
                 superheat=20, eta_pump=0.8):
                     
        vars = """Q_condenser,W
Q_boiler,W
Q_evaporator,W
W_pump,W
COP,W/W
m_dot,kg/s
m_refrig,kg/s""".split()
        self.vars = [var.split(',')[0] for var in vars]
        self.units = [var.split(',')[1] for var in vars]
        for var in self.vars:
            self.__setattr__(var,0)
            
        self.T_hot = T_hot
        self.T_reject = T_reject
        self.T_cold = T_cold
        self.m_dot = m_dot
        self.eta_nozzle = eta_nozzle
        self.eta_mixer = eta_mixer
        self.eta_diffuser = eta_diffuser
        self.superheat = superheat
        self.eta_pump = eta_pump
    
        self.labels=range(1,12)+["1'","5'","7'"]
        self.points={}
        for i in self.labels:
            self.points[i] = CP.AbstractState("HEOS","Water")
        
        
        
    def update(self,debug=False):
        # Condenser
        self.points[9].update(CP.QT_INPUTS,0,C2K(self.T_reject))
        self.points[8].update(CP.PQ_INPUTS,self.points[9].p(),1)
        
        # Boiler
        self.points[3].update(CP.QT_INPUTS,1,C2K(self.T_hot))
        self.points[4].update(CP.PT_INPUTS,
            self.points[3].p(),self.points[3].T()+self.superheat)
        self.points[2].update(CP.PQ_INPUTS,
            self.points[3].p(),0)
        # Pump
        p1prime = CP.AbstractState("HEOS","Water")
        p1prime.update(CP.PSmass_INPUTS,
                       self.points[2].p(),self.points[9].smass())
        self.points["1'"] = p1prime
        dhprime = p1prime.hmass() - self.points[9].hmass()
        dh = dhprime / self.eta_pump
        self.points[1].update(CP.HmassP_INPUTS,
            self.points[9].hmass() + dh, self.points[2].p())
        
        # Evaporator and expansion valve
        self.points[11].update(CP.QT_INPUTS,1,C2K(self.T_cold))
        self.points[10].update(CP.HmassP_INPUTS,
            self.points[9].hmass(),self.points[11].p())
        
        # Now the fun begins!
        # Ideal, isentropic expansion through nozzle
        p5prime = CP.AbstractState("HEOS","Water")
        p5prime.update(CP.PSmass_INPUTS,
                       self.points[11].p(),self.points[4].smass())
        self.points["5'"] = p5prime
        dhprime = p5prime.hmass() - self.points[4].hmass()
        dh = self.eta_nozzle * dhprime
        KE5 = -dh
        self.points[5].update(CP.HmassP_INPUTS,
            self.points[4].hmass() + dh, self.points[11].p())
        
        # Entrainment
        p7prime = CP.AbstractState("HEOS","Water")
        h_6 = self.points[5].hmass()
        self.points[6].update(CP.HmassP_INPUTS,
                h_6, self.points[5].p())
        
        if debug:
            print "h4 = ", self.points[4].hmass()
            print "KE5 = ", KE5
            print "guess", h_6
        
        for i in range(11):
            p7prime.update(CP.PSmass_INPUTS,
                self.points[8].p(), self.points[6].smass())
            
            # Equation 4
            dhprime = p7prime.hmass() - self.points[6].hmass()
            dh = dhprime / self.eta_diffuser
            self.points[7].update(CP.HmassP_INPUTS,
                self.points[6].hmass() + dh, self.points[8].p())
            # Equation 2
            KE6 = self.points[7].hmass() - h_6
            # Equation 3
            m_out = self.eta_mixer * self.m_dot * KE5 / KE6
            self.m_refrig = m_out - self.m_dot
            # Equation 1
            h_6 = (self.m_dot * self.points[4].hmass() \
                + self.m_refrig * self.points[11].hmass()) / m_out - KE6
            self.points[6].update(CP.HmassP_INPUTS,
                h_6, self.points[5].p())
            if debug:
                print i, h_6, self.points[7].hmass(), KE6, self.m_refrig, m_out
        self.points["7'"] = p7prime
        if debug: print
        # Heat and work flows
        self.Q_condenser = m_out * (self.points[9].hmass() - self.points[7].hmass())
        self.Q_boiler = self.m_dot * (self.points[4].hmass() - self.points[1].hmass())
        self.Q_evaporator = self.m_refrig * (self.points[11].hmass() - self.points[10].hmass())
        self.COP=self.Q_evaporator/self.Q_boiler
        self.W_pump = self.m_dot * (self.points[1].hmass() - self.points[9].hmass())
        
    def __repr__(self):
        keys = [CP.iT, CP.iP, CP.iHmass, CP.iQ, CP.iSmass]
        states = []
        for k,p in self.points.iteritems():
            s = [safe_keyed_output(p,i) for i in keys]
            states.append([k,]+s)
        statestr=tabulate.tabulate(states,["i","T","P","Hmass","Q","Smass"])
        data=[]
        for var,unit in zip(self.vars,self.units):
            data.append((var,unit,self.__getattribute__(var)))
        varstr = tabulate.tabulate(data)
        return statestr+"\n\n"+varstr
        
    
    def display(self):
        import matplotlib.pyplot as plt
        for stream,ii in [("refrig", [9,10,11]),
                          ("mix", [7,8,9]), #6, ...
                          ("motive", [9,1,2,3,4]), #... ,5
                          ("pump", [9,"1'"]),]:
                          #("nozzle", [4,"5'"]),
                          #("diffuser", [6, "7'"])]:
            hh,TT = zip(*[(self.points[i].hmass(),self.points[i].T()) for i in ii])
            plt.plot(hh,TT)
        state = CP.AbstractState("HEOS","Water")
        # Fix nozzle stream
        for i1,i2,c in [(4,5,'r'),(4,"5'",'r--'),(6,7,'g'),(6,"7'",'g--')]:
            p1,p2 = self.points[i1],self.points[i2]
            ds = p2.smass() - p1.smass()
            dp = p2.p() - p1.p()
            HH = []
            TT = []
            for x in np.linspace(0,1):
                s = ds * x + p1.smass()
                p = dp * x + p1.p()
                state.update(CP.PSmass_INPUTS,p,s)
                HH.append(state.hmass())
                TT.append(state.T())
            plt.plot(HH,TT,c)
        
        # Saturation curves
        #state.update(CP.QT_INPUTS,0,state.Tmin())
        #pmin = state.p()
        #pmax = state.p_critical()
        tmin = state.Tmin()
        tmax = state.T_critical()
        lh,lT = [],[]
        vh,vT = [],[]
        #for P in np.linspace(pmin, pmax):
        for T in np.linspace(tmin,tmax):
            #state.update(CP.PQ_INPUTS, P, 0)
            state.update(CP.QT_INPUTS, 0, T)
            lh.append(state.hmass())
            lT.append(state.T())
            #state.update(CP.PQ_INPUTS, P, 1)
            state.update(CP.QT_INPUTS, 1, T)
            vh.append(state.hmass())
            vT.append(state.T())
        plt.plot(lh,lT,'--')
        plt.plot(vh,vT,'--')
        
    def plotInternalCurves(self):
        import matplotlib.pyplot as plt
        plt.figure()
        for c,Q in [(self.getEvaporatorCurve,self.Q_evaporator),
                    (self.getCondenserCurve,self.Q_condenser),
                    (self.getBoilerCurve,self.Q_boiler)]:
            curve = c()
            q = np.linspace(0,Q*1.05)
            T = curve.T(q)            
            plt.plot(q,T)            
        
    def getBoilerCurve(self):
        return HRHX_integral_model.waterStream(self.points[3].p(),
                                               self.points[1].hmass(),
                                               self.m_dot)
    def getCondenserCurve(self):
        return HRHX_integral_model.waterStream(self.points[9].p(),
                                               self.points[7].hmass(),
                                               self.m_dot + self.m_refrig)
    def getEvaporatorCurve(self):
        return HRHX_integral_model.waterStream(self.points[11].p(),
                                               self.points[10].hmass(),
                                               self.m_refrig)

class EjectorSystem(object):
    def __init__(self):
        self.boiler = HRHX_integral_model.streamExample1(T_inlet=250,mdot=11,cp=4215)
        self.condenser = HRHX_integral_model.streamExample1(45,200,4215)
        self.evaporator = HRHX_integral_model.streamExample1(12,3,4215)
    def calcUA(self,ec,plot=False):
        hrhx1 = HRHX_integral_model.counterflow_integrator(
            ec.getBoilerCurve(), self.boiler)
        DeltaT1,epsilon1,UA1 = hrhx1.calcUA2(ec.Q_boiler)
        if plot:
            HRHX_integral_model.plotFlow(hrhx1,None,ec.Q_boiler)
            plt.xlabel("Heat (W)")
            plt.ylabel("Temperature ($^\circ$C)")
            plt.title("Boiler $\Delta T$ = {:g}, UA = {:g}".format(DeltaT1,UA1))
            
        hrhx2 = HRHX_integral_model.counterflow_integrator(
            self.condenser, ec.getCondenserCurve())
        DeltaT2,epsilon2,UA2 = hrhx2.calcUA2(-ec.Q_condenser)
        if plot:
            HRHX_integral_model.plotFlow(hrhx2,None,-ec.Q_condenser)
            plt.title("Condenser $\Delta T$ = {:g}, UA = {:g}".format(DeltaT2,UA2))
            
        hrhx3 = HRHX_integral_model.counterflow_integrator(
            ec.getEvaporatorCurve(), self.evaporator)
        DeltaT3,epsilon3,UA3 = hrhx3.calcUA2(ec.Q_evaporator)
        if plot:
            HRHX_integral_model.plotFlow(hrhx3,None,ec.Q_evaporator)
            plt.title("Evaporator $\Delta T$ = {:g}, UA = {:g}".format(DeltaT3,UA3))
        
        return (UA1 + UA2 + UA3),[DeltaT1,DeltaT2,DeltaT3]

if __name__ == "__main__":
    x0 = (200,60,5,1)
    ec = EjectorCycle(*x0)
    ec.update()
    print ec
    ec.display()
    ec.plotInternalCurves()
    plt.show()
    
    es = EjectorSystem()
    UA0 = es.calcUA(ec)[0]
    plt.show()
    
    def objective(x):
        print "Objective at x = ", x
        try:
            ec = EjectorCycle(*x)
            ec.update()
            UA, DTs = es.calcUA(ec)
        except:
                return 0
        if np.less(DTs,0).any():
            return 0
        else:
            return -ec.Q_evaporator
    def constraint0(x):
        try:
            ec = EjectorCycle(*x)
            ec.update()
            return UA0 - es.calcUA(ec)[0]
        except:
            return -1
    def constraint1(x):
        return 400 - x[0]
    def constraint2(x):
        return 400 - x[1]
    def constraint3(x):
        return 400 - x[2]
    def constraint4(x):
        return 5 - x[3]
    def constraint5(x):
        return x[0]
    def constraint6(x):
        return x[1]
    def constraint7(x):
        return x[2]
    def constraint8(x):
        return x[3]
    constraints = [{"type":"ineq","fun":c} for c in [constraint0,constraint1,constraint2,constraint3,constraint4,constraint5,constraint6,constraint7,constraint8]]
    opt = scipy.optimize.minimize(objective,x0,constraints=constraints)