# -*- coding: utf-8 -*-
"""
Created on Wed May 25 19:54:52 2016

@author: nfette
"""
import scipy.integrate
import scipy.optimize
import matplotlib.pyplot as plt
import numpy as np

class stream(object):
    def setQ(self,Q):
        pass
    def q(self,T):
        pass
    def T(self,q):
        pass

class streamExample1(object):
    def __init__(self, T_inlet=0, mdot=1.0, cp=1.0):
        self.C = mdot * cp
        self.T_inlet = T_inlet
    def q(self,T):
        return self.C * (T - self.T_inlet)
    def T(self,q):
        return self.T_inlet + q / self.C

class streamExample2(object):
    """Pure fluid boiling at constant temperature.
    Input the inlet enthalpy relative to h = 0 at saturated liquid."""
    def __init__(self,h_inlet,T_sat,mdot,cp,hsatv):
        self.h_inlet = h_inlet
        self.T_sat = T_sat
        self.mdot = mdot
        self.cp = cp
        self.hsatl = 0
        self.hsatv = hsatv
        self.T = np.vectorize(self._T)
        self.q = np.vectorize(self._q)
    def _q(self,T):
        sensible = (T - self.T_sat) * self.cp
        test = (self.mdot > 0) * (T > self.T_sat) + (self.mdot < 0) *  (T >= self.T_sat)
        latent = self.hsatv if test else self.hsatl
        dh = sensible + latent - self.h_inlet
        return dh * self.mdot
    def _T(self,q):
        dh = q / self.mdot
        h = dh + self.h_inlet
        if h < self.hsatl:
            # h = hsatl + cp * (T - T_sat)
            return self.T_sat + (h - self.hsatl) / self.cp
        elif h > self.hsatv:
            # h = hsatl + cp * (T - T_sat)
            return self.T_sat + (h - self.hsatv) / self.cp
        else:
            return self.T_sat
    
        
class counterflow_integrator(object):
    def __init__(self, cold, hot, useHotT=False):
        self.cold = cold
        self.hot = hot
        self.useHotT = useHotT
        self.calcQmax()
    def calcUA(self,Q):
        func = lambda(q):1./(self.hot.T(Q-q)-self.cold.T(q))
        return scipy.integrate.quad(func,0,Q)[0]
    def calcQ(self,UA):
        func = lambda(Q):(self.calcUA(Q)-UA)**2
        constraints = [{"type":"ineq",
                       "fun":lambda(Q):Q},
                        {"type":"ineq",
                       "fun":lambda(Q):self.Qmax-Q}]
        return scipy.optimize.minimize(func,0,constraints=constraints).x[0]
    def calcQmax(self,extra=False,brute=True):
        if self.useHotT:
            func = lambda(Q):Q+self.cold.q(self.hot.T(Q))
        else:
            func = lambda(Q):Q+self.hot.q(self.cold.T(Q))
        qc = min(self.cold.q(self.hot.T(0)), self.hot.q(self.cold.T(0)))
        print "qc = ",qc
        lim = lambda(Q):qc-Q
        constraints = [{"type":"ineq",
                       "fun":lambda(Q):Q,
                       "jac":lambda(Q):1},
                       {"type":"ineq",
                        "fun":lim,
                        "jac":lambda(Q):-1}]
        if brute:
            opt = scipy.optimize.differential_evolution(func,[(0,qc)])
            print opt
        else:
            opt = scipy.optimize.minimize(func,0,constraints=constraints)
        self.Qmax = np.asscalar(opt.fun)
        if extra:
            return opt
        else:
            return self.Qmax

def plotFlow(ci):
    QQ = np.linspace(0,ci.Qmax)
    Tcold = ci.cold.T(QQ)
    Thot1 = ci.hot.T(QQ)
    plt.figure()
    plt.plot(QQ,Tcold,'b.',label="cold")
    for Q in [0., ci.Qmax * 0.5,ci.Qmax * 0.99]:
        UA=ci.calcUA(Q)
        #plt.plot(QQ,Tcold,'b')
        plt.plot(Q-QQ,Thot1,'.',label="Q={:.3f},UA={:.3f}".format(Q,UA))
    #plt.legend(loc='best')
    
def plotNN(ci):    
    NNs = np.arange(1,10)
    QQ = ci.Qmax * (1 - np.exp(-NNs))
    UA = map(ci.calcUA,QQ)
    plt.figure()
    plt.plot(UA,QQ,'ko-')
    plt.xlabel("UA")
    plt.ylabel("Q")

if __name__ == "__main__":
    cold = streamExample1(0)
    hot = streamExample1(1,-1)
    ci = counterflow_integrator(cold, hot)
    Qmax = ci.Qmax
    print "Qmax=",Qmax
    UA=ci.calcUA(Qmax/2)
    print "UA=",UA
    
    plotFlow(ci)
    plotNN(ci)
    
    c2 = streamExample2(-5.,100.,1.,1.,10.)
    h2 = streamExample1(120.,-1.5)
    ci2 = counterflow_integrator(c2,h2)
    print ci2.Qmax
    
    plotFlow(ci2)
    plotNN(ci2)
    plt.show()
    
    c3 = streamExample2(-5.,100.,1.,1.,10.)
    h3 = streamExample2(15.,100.,-1.,1.,10.)
    ci3 = counterflow_integrator(c3,h3)
    print ci3.Qmax
    
    plotFlow(ci3)
    plotNN(ci3)
    plt.show()
    
    QQ=np.arange(0,35)
    Tc=c3.T(QQ)
    Th=h3.T(QQ)
    TT=np.arange(95,110)
    Qc=c3.q(TT)
    Qh=h3.q(TT)
    plt.figure();plt.plot(QQ,Tc,QQ,Th);plt.plot(Qc,TT,'.',Qh,TT,'.')