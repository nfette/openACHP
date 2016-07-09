# -*- coding: utf-8 -*-
"""
Created on Wed May 25 19:54:52 2016

@author: nfette
"""
import scipy.integrate
import scipy.optimize
import matplotlib.pyplot as plt
import numpy as np
from hw2_1 import CelsiusToKelvin as C2K, KelvinToCelsius as K2C
from exceptions import ValueError

class stream(object):
    def setQ(self,Q):
        pass
    def q(self,T):
        pass
    def T(self,q):
        pass

class streamExample1(object):
    def __init__(self, T_inlet=0, mdot=1.0, cp=1.0):
        self.mdot = mdot
        self.cp = cp
        self.T_inlet = T_inlet
    def q(self,T):
        C = self.mdot * self.cp
        return C * (T - self.T_inlet)
    def T(self,q):
        C = self.mdot * self.cp
        return self.T_inlet + q / C

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
        if T < self.T_sat:
            latent = self.hsatl
            h_out = latent + sensible
        elif T > self.T_sat:
            latent = self.hsatv
            h_out = latent + sensible
        else:
            if self.h_inlet < self.hsatl:
                h_out = self.hsatl
            elif self.h_inlet > self.hsatv:
                h_out = self.hsatv
            else:
                h_out = self.h_inlet
        
        dh = h_out - self.h_inlet
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

class waterStream(object):
    """A class for phase change of a given fluid. The specific heats of
    saturated water liquid and vapor at equilibrium are different,
    so the streamExample2 is insufficient to represent this case."""
    def __init__(self,P,h_in,mdot,fluid='water'):
        import CoolProp
        
        self.P = P
        self.h_in = h_in
        self.mdot = mdot
        self.fluid = fluid
        f = CoolProp.AbstractState('HEOS',fluid)
        f.update(CoolProp.PQ_INPUTS,self.P,0)
        h_liq = f.hmass()
        f.update(CoolProp.PQ_INPUTS,self.P,1)
        h_vap = f.hmass()
        # Determine what enthalpy to give at saturation temperature, since
        # at saturation temperature, rounding depends whether
        # the user intends heating or cooling.
        if h_in < h_liq:
            self.h_sat = h_liq
        elif h_in > h_vap:
            self.h_sat = h_vap
        else:
            self.h_sat = h_in
        
        #f.update(CoolProp.QT_INPUTS,0,f.Ttriple())
        f.update(CoolProp.PT_INPUTS,P,f.Tmin())
        h_min = f.hmass()
        f.update(CoolProp.PT_INPUTS,P,f.Tmax())
        h_max = f.hmass()
        qlim1,qlim2 = self.mdot * (h_max - self.h_in), self.mdot * (h_min - self.h_in)
        self.qmin = min(qlim1,qlim2)
        self.qmax = max(qlim1,qlim2)
        self.Tmin = K2C(f.Tmin())
        self.Tmax = K2C(f.Tmax())
        self.q = np.vectorize(self._q)
        self.T = np.vectorize(self._T)
            
    def _q(self,T):
        from CoolProp.CoolProp import PropsSI
        if T < self.Tmin:
            return 0
        elif T > self.Tmax:
            return self._q(self.Tmax)
        else:
            try:
                # CoolProp raises an exception if this is the saturation temp.
                #h_out = CP.PropsSI('H','P',self.P,'T',C2K(T),self.fluid)
                h_out = PropsSI('H','P',self.P,'T',C2K(T),self.fluid)
            except:
                h_out = self.h_sat
        return self.mdot * (h_out - self.h_in)

    def _T(self,q):        
        from CoolProp.CoolProp import PropsSI
        if q < self.qmin:
            q = self.qmin
        elif q > self.qmax:
            q = self.qmax
        
        h_out = q / self.mdot + self.h_in
        return K2C(PropsSI('T','P',self.P,'H',h_out,self.fluid))

class counterflow_integrator(object):
    """Change in progress:
    In theoretical document, sign convention was that q > 0 for both cold
    and hot streams. However, discontinuities in q(T) for phase change turned
    out to be evil, and so it will be simpler to adopt a different convention:
    q > 0: heat is transferred into stream (cold stream is heating up)
    q < 0: heat is transferred out of stream (hot stream is cooling down)
    """
    def __init__(self, cold, hot, useHotT=False):
        self.cold = cold
        self.hot = hot
        self.useHotT = useHotT
        # Functions for Qmax
        self.func1 = lambda(Q):Q+self.cold.q(self.hot.T(-Q))
        self.func2 = lambda(Q):Q-self.hot.q(self.cold.T(Q))
        
        self.calcQmax()
    def calcUA(self,Q,eff=False):
        # OLD
        #func = lambda(q):1./(self.hot.T(Q-q)-self.cold.T(q))
        # Q > is total heat transferred into cold stream, and q is local cum.
        func = lambda(q):1./(self.hot.T(q-Q)-self.cold.T(q))
        ua = scipy.integrate.quad(func,0,Q)[0]
        epsilon = Q / self.Qmax
        if epsilon > 1:
            raise ValueError("Q given [{}] is higher than Q maximum [{}];"\
            " effectiveness [{}] > 1.".format(Q,self.Qmax,epsilon))
        if eff:            
            return ua, epsilon
        else:
            return ua
    def calcQ(self,UA):
        func = lambda(Q):(self.calcUA(Q)-UA)**2
        constraints = [{"type":"ineq",
                       "fun":lambda(Q):Q},
                        {"type":"ineq",
                       "fun":lambda(Q):self.Qmax-Q}]
        return scipy.optimize.minimize(func,0,constraints=constraints).x[0]
    def calcQmax(self,extra=False,brute=False):
        # Preliminary Max Q based on inlet temperatures only
        qc = min(self.func1(0),self.func2(0))
        #print "qc = ",qc
        lim = lambda(Q):qc-Q
        constraints = [{"type":"ineq",
                       "fun":lambda(Q):Q,
                       "jac":lambda(Q):1},
                       {"type":"ineq",
                        "fun":lim,
                        "jac":lambda(Q):-1}]
        if brute:
            opt1 = scipy.optimize.differential_evolution(self.func1,[(0,qc)])
            opt2 = scipy.optimize.differential_evolution(self.func2,[(0,qc)])
            #print opt1
            #print opt2
        else:
            opt1 = scipy.optimize.minimize(self.func1,0,constraints=constraints)
            opt2 = scipy.optimize.minimize(self.func2,0,constraints=constraints)
            #print opt1
            #print opt2
        self.Qmax = min(np.asscalar(opt1.fun),np.asscalar(opt2.fun))
        #self.Qmax = np.asscalar(opt.x)
        if extra:
            return opt1
        else:
            return self.Qmax

def plotFlow(ci,figure=None,Qactual=None):
    QQ = np.linspace(0,ci.Qmax)
    Tcold = ci.cold.T(QQ)
    Thot1 = ci.hot.T(-QQ)
    if not figure:
        plt.figure()
    plt.plot(QQ,Tcold,'b.-',label="cold")
    for Q,c in zip([0., ci.Qmax * 0.5,ci.Qmax * 0.99,Qactual],['red','orange','yellow','black']):
        if Q is not None:
            #UA=ci.calcUA(Q)
            UA = 0
            #plt.plot(QQ,Tcold,'b')
            #plt.plot(Q-QQ,Thot1,'.-',label="Q={:.3f},UA={:.3f}".format(Q,UA),color=c)
            plt.plot(Q-QQ,Thot1,'.-',color=c)
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
    if True:
        cold = streamExample1(0)
        hot = streamExample1(1,1)
        ci = counterflow_integrator(cold, hot)
        Qmax = ci.Qmax
        print "Qmax=",Qmax
        UA=ci.calcUA(Qmax/2)
        print "UA=",UA
        
        plotFlow(ci)
        plotNN(ci)
        
        c2 = streamExample2(-5.,100.,1.,1.,10.)
        h2 = streamExample1(120.,1.5)
        ci2 = counterflow_integrator(c2,h2)
        print ci2.Qmax
        
        plotFlow(ci2)
        plotNN(ci2)
        plt.show()
        
        c3 = streamExample2(-5.,100.,1.,1.,10.)
        h3 = streamExample2(15.,100.,1.,1.,10.)
        # The first result with local optimization routine is wrong,
        # because the distance function has discontinuities.
        ci3 = counterflow_integrator(c3,h3)
        # However, a global optimization finds the correct basin.
        ci3.calcQmax(extra=False,brute=True)
        print ci3.Qmax
        
        plotFlow(ci3)
        plotNN(ci3)
        plt.show()
        
    
    water1=waterStream(101325,100,1)
    TT=np.linspace(20,120)
    QQ=water1.q(TT)
    Q2=np.linspace(min(QQ),max(QQ))
    T2=water1.T(Q2)
    plt.figure()    
    plt.plot(QQ,TT,'.',Q2,T2,'-')
    water2=waterStream(15e5,2.8e6,1)
    ci4 = counterflow_integrator(water1,water2,useHotT=True)
    plotFlow(ci4)
    plt.show()
    plotNN(ci4)
    plt.show()