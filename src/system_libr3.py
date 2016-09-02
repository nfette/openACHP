# -*- coding: utf-8 -*-
"""
Created on Sat Apr 30 03:26:36 2016

This file solves a thermal system with a chiller.

@author: nfette
"""

import numpy as np
import CoolProp.CoolProp as CP
import tabulate
import scipy.optimize
from scipy.special import expit, logit

from hw2_1 import CelsiusToKelvin as C2K
from hw2_1 import KelvinToCelsius as K2C

import libr3
from HRHX_integral_model import streamExample1 as se1, \
    counterflow_integrator as hxci, \
    plotFlow as plotFlow


def makeSystem():
    return system1(
        heatStream=se1(120,0.3,4179),
        absorberRejectStream=se1(15,1,4179),
        condRejectStream=se1(35,1,4179),
        coldStream=se1(12,1,4179))

def makeChiller(x):
    m_pump,T_evap,T_cond,x1,x2 = x
    chiller = libr3.ChillerLiBr1(T_evap, T_cond, x1, x2, m_pump=m_pump)
    chiller.iterate1()
    return chiller

def boundmap(xmin,xmax,x):
    alpha = (x - xmin) / (xmax - xmin)
    return logit(alpha)
def unmap(xmin,xmax,z):
    alpha = expit(z)
    return xmin + (xmax - xmin) * alpha
    
class varBounds(object):
    """Provides a map for re-parameterizing bounded variables.
    Bounds are bounds on domain X.
    Forward maps x in X to y in Y = R^n.
    Backward maps y in Y to x in X."""
    def __init__(self):
        self.mmin = 0.0
        self.mmax = 0.3
        self.Temin = 1.0
        self.Temax = 6.0
        self.Tcmin = 10.
        self.Tcmax = 65.
        self.xmin = 0.4
        self.xmax = 0.7
    def forward(self,x):
        a = boundmap(self.mmin,self.mmax,x[0])
        b = boundmap(self.Temin,self.Temax,x[1])
        c = boundmap(self.Tcmin,self.Tcmax,x[2])
        d = boundmap(self.xmin,self.xmax,x[3])
        e = boundmap(x[3],self.xmax,x[4])
        return [a,b,c,d,e]
    def backward(self,y):
        a = unmap(self.mmin,self.mmax,y[0])
        b = unmap(self.Temin,self.Temax,y[1])
        c = unmap(self.Tcmin,self.Tcmax,y[2])
        d = unmap(self.xmin,self.xmax,y[3])
        e = unmap(d,self.xmax,y[4])
        return [a,b,c,d,e]
        
class system1:
    def __init__(self,heatStream,absorberRejectStream,condRejectStream,coldStream):
        self.heatStream = heatStream
        self.absorberRejectStream = absorberRejectStream
        self.condRejectStream = condRejectStream
        self.coldStream = coldStream
        
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
        print("Objective: x = {}".format(x))
        mdot,T_evap,T_cond,x1,x2 = x
        try:
            UA=self.calcUA(x)
            Q = self.chiller.Q_evap_heat
            print(" -> Q = {}, UA = {}".format(Q,UA))
            return -Q
        except:
            print("An error occurred in calcUA(). Reporting zeros.")
            return 0
                  
    def constraint1(self,x,UAgoal):
        # Result is zero iff: UA is the given value.
        print("Constraint1: x={}, UAgoal={}".format(x,UAgoal))
        try:
            UA = self.calcUA(x)
            print("UA at x = {}".format(UA))
            return UAgoal - UA
        except:
            print("An error occurred in calcUA(). Reporting zeros.")
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
                     dict(type='ineq',fun=self.constraint2,jac=lambda x:[0,0,0,-1,1]),
                     dict(type='ineq',fun=self.constraint3,jac=lambda x:[0,-1,1,0,0]),
                     dict(type='ineq',fun=self.constraint4,jac=lambda x:[1,0,0,0,0]),
                     dict(type='ineq',fun=self.constraint5,jac=lambda x:[0,0,0,1,0]),
                     dict(type='ineq',fun=self.constraint6,jac=lambda x:[0,0,0,0,-1]),
                     dict(type='ineq',fun=self.constraint7,jac=lambda x:[0,1,0,0,0]),
                     dict(type='ineq',fun=self.constraint8,jac=lambda x:[0,0,-1,0,0]),
                     dict(type='ineq',fun=self.constraint9,jac=lambda x:[-1,0,0,0,0])
                    ]
                
        opt = scipy.optimize.minimize(self.objective,x0,constraints=constraints,
                                      method='COBYLA',
                                      options=dict(disp=True,maxiter=10))
        return opt

class pair(object):
    def __init__(self,sys,chiller):
        #self.coldStream = self.makeColdStream(chiller.Q_evap_heat)
        #self.absorberRejectStream = self.makeRejectStream(32, 5, chiller.Q_abs_total)
        #self.condRejectStream = self.makeRejectStream(20, 5, chiller.Q_condenser_reject)
        self.sys = sys
        self.chiller = chiller
        
        self.hxs = []
        self.hxs.append(("gen",
            hxci(chiller.getGeneratorStream(),sys.heatStream,useHotT=True),
            chiller.Q_gen_total))
        self.hxs.append(("abs",
            hxci(sys.absorberRejectStream,chiller.getAbsorberStream()),
            chiller.Q_abs_total))
        self.hxs.append(("cond",
            hxci(sys.condRejectStream,chiller.getCondenserStream()),
            chiller.Q_condenser_reject))
        self.hxs.append(("evap",
            hxci(chiller.getEvaporatorStream(),sys.coldStream),
            chiller.Q_evap_heat))
        
        self.data = []
        self.totalUA = 0
        for name, ci, Q in self.hxs:
            ci.calcQmax()
            deltaT, epsilon, UA = ci.calcUA2(Q)
            self.data.append((name, deltaT, epsilon, UA, Q))
            self.totalUA += UA
    
    def __repr__(self):
        result = tabulate.tabulate(self.data,"name deltaT epsilon UA Q".split())
        result += "\ntotalUA = {}".format(self.totalUA)
        return result
    
    def display(self):
        import matplotlib.pyplot as plt
        figs = []
        for name, ci, Q in self.hxs:
            fig=plotFlow(ci, None, Q)
            figs.append(fig)
            fig.gca().set_title(name)
            plt.xlabel("Heat (W)")
            plt.ylabel("Temperature (deg C)")
        return figs
    
    def constraints(self,UAgoal):
        c = []
        
        # Result is >=0 iff: x2 > x1
        #c.append((chiller.x2 - chiller.x1,0.01))
        
        # Result is >=0 iff: T_cond > T_evap
        #c.append((chiller.T_cond - chiller.T_evap,1))

        # Result is >=0 iff: mdot > 0
        #c.append((chiller.m_pump,0.01))
        
        # Result is >=0 iff: x1 > 0
        #c.append((chiller.x1,0.01))
        
        # Result is >=0 iff: x2 < 1
        #c.append((0.7 - chiller.x2,0.01))
    
        # Result is >=0 iff: T_evap > 0
        #c.append((chiller.T_evap,0.1))
        
        # Result is >=0 iff: T_evap > 0
        #c.append((65 - chiller.T_cond,0.1))
        
        # Result is >=0 iff: T_evap > 0
        #c.append((0.3 - chiller.m_pump,0.01))

        for name,deltaT,epsilon,UA,Q in self.data:
            c.append((deltaT,0.1))
            #c.append(deltaT)
            
        # Result is zero iff: UA is the given value.
        if self.totalUA is not np.inf:
            c.append((UAgoal - self.totalUA, 0.01 * UAgoal))
        #c.append(UAgoal - self.totalUA)

        return c
    
    def objective(self,UAgoal):
        result = -self.chiller.Q_evap_heat
        for c,mu in self.constraints(UAgoal):
            result *= expit(c/mu)
        return result

class objectivePair(object):
    def __init__(self,sys,UAgoal):
        self.sys = sys
        self.UAgoal = UAgoal
        self.x=[]
        self.bounds = varBounds()
        self.plot = progressPlot()
    def __call__(self,y):
        x = self.bounds.backward(y)
        print("objective at x = ", x)
        self.x.append(x)
        
        try:
            chiller = makeChiller(x)
            self.plot.updateQ(chiller.Q_evap_heat)
            p = pair(self.sys, chiller)
            self.plot.updateUA([row[3] for row in p.data])
            result = p.objective(self.UAgoal)
            print("Returning ", result)
            return result
        except Exception as e:
            print(e)
            print("Returning zero anyway")
            return 0

def main():
    import matplotlib.pyplot as plt
    import numpy as np
    
    sys = makeSystem()
    
    #x0 = [m_pump, T_evap, T_cond, x1, x2]
    UAgoal = 10000
    x0 = [0.067814426902490302, 2.035419430770669, 52.071417867824564, 0.51537879930261832, 0.61044406504494919]
    x0 = [0.067912042024285593, 1.0000000420821966, 52.061185313461571, 0.51227487505032543, 0.6093890922208236]
    #UAgoal = 1000
    #x0 = [0.029239052322855189, 2.5596320340218628, 52.383189871781816, 0.57568178054004149, 0.61712740762371499]
    #x0 = [0.019500209749505271, 1.0000000026568656, 53.095537722738428, 0.55188715550421308, 0.59664237804180931]
    #x0 = [0.012743820293376019, 1.0000000026568656, 53.01139267358721, 0.54142012257312411, 0.6125033847474235]
    #x0 = [0.012743828571511037, 1.0000000026590252, 53.008775297974509, 0.54142083075023628, 0.61250377548300572]
    #UAgoal = 100
    #x0 = [0.001978486398852077, 1.8724708108120764, 49.776208068337269, 0.57164104762999846, 0.61708756538489995]
    #x0 = [0.0013539934436546737, 1.0000050517501358, 53.173634288614451, 0.54343850911212377, 0.61262006880404485]
    #x0 = [0.0013540160549637599, 1.0000029445505649, 53.171993101646194, 0.54343663586360358, 0.61261476185248154]
    #x0 = [0.0013540064551193383, 1.0000029445505414, 53.169536604904081, 0.54345636471281844, 0.61262581255195125]
    # 10
    #x0 = [0.00012243761622050806, 1.5744486090660881, 53.296690938447185, 0.53538141034516507, 0.61112929898671542]


    chiller = makeChiller(x0)
    
    p = pair(sys,chiller)
    print(p)
    print(p.hxs[0][1].cold)
    figs = p.display()

    if False:
        fig = plt.figure()
        plt.title("Generator")
        for mdot in np.linspace(0.01,0.5,10):
            try:
                x1 = x0[:]
                x1[0] = mdot
                ch1 = makeChiller(x1)
                p1 = pair(sys,ch1)
                print(p1)
                plotFlow(p1.hxs[0][1],fig,p1.hxs[0][2])
            except Exception as e:
                print(e)
                break
    plt.show()
    
    
    op = objectivePair(sys,UAgoal)
    y0 = op.bounds.forward(x0)
    print(y0)
    try:
        opt = scipy.optimize.minimize(op,y0,
            method='Powell',options=dict(maxftol=0.001,maxfev=100))
        #options=dict(maxiter=10),tol=0.1
    except Exception as e:
        print(e)
        opt = None
        
    return sys, chiller, p, figs, op, opt

class progressPlot(object):
    def __init__(self):
        #http://stackoverflow.com/questions/10944621/dynamically-updating-plot-in-matplotlib
        import matplotlib.pyplot as plt
        # For Q_evap progress
        self.fig1,self.ax1 = plt.subplots()
        self.line1, = self.ax1.plot([0,1],[0,1])
        self.Q = []
        
        # For UA progress
        self.fig2,self.ax2 = plt.subplots()    
        self.lines = [self.ax2.plot([0,1],[0,1])[0] for i in range(4)]
        self.UAs=[[],[],[],[]]
        
    def updateQ(self,q):
        self.Q.append(q)
        self.line1.set_data(range(len(self.Q)),self.Q)
        self.ax1.relim()
        self.ax1.autoscale_view()
        self.fig1.canvas.draw()
        self.fig1.canvas.flush_events()
        
    def updateUA(self,row):
        for j,ua in enumerate(row):
            self.UAs[j].append(row[j])
            self.lines[j].set_data(range(len(self.UAs[j])),self.UAs[j])
        self.ax2.relim()
        self.ax2.autoscale_view()
        self.fig2.canvas.draw()
        self.fig2.canvas.flush_events()
    
    def postProcess(self,sys,op):
        """Input a system1 and objectivePair.
        Then loop through the x values in op and plot the UA and Q."""
        for i,x in enumerate(op.x):
            ch=makeChiller(x)
            self.updateQ(ch.Q_evap_heat)
            p=pair(sys,ch)
            self.updateUA([row[3] for row in p.data])

if __name__ == "__main__":
    sys, chiller, p, figs, op, opt = main()
    print(opt)
    xopt = op.bounds.backward(opt.x)
    print(xopt)
    ch1 = makeChiller(xopt)
    p1 = pair(sys,ch1)
    print(p1)
    p1.display()
    