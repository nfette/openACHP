# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 15:48:09 2016

@author: nfette
"""
import tabulate
import HRHX_integral_model
import numpy as np
import matplotlib.pyplot as plt

class CarnotEngine(object):
    """Temperature in K"""
    def __init__(self,T_hot,T_cold,m_dot=1.,h_fg=100.):
        self.T_hot = T_hot
        self.T_cold = T_cold
        self.m_dot = m_dot
        self.h_fg = h_fg
        self.Q_hot = m_dot * h_fg
        self.eta = 1. - T_cold / T_hot
        self.W_net = - self.eta * self.Q_hot
        self.Q_cold = - self.Q_hot - self.W_net
    def __repr__(self):
        names = """T_hot T_cold m_dot h_fg Q_hot eta W_net Q_cold""".split()
        result = []
        for i in names:
            result.append(self.__getattribute__(i))
        return tabulate.tabulate(zip(names,result))
    def getCondenserStream(self):
        deltaH = - self.Q_cold / self.Q_hot * self.h_fg
        return HRHX_integral_model.streamExample2(
            deltaH,self.T_cold,self.m_dot,1.0,deltaH)
    def getBoilerStream(self):
        return HRHX_integral_model.streamExample2(
            0,self.T_hot,self.m_dot,1.0,self.h_fg)

class CarnotCooler(object):
    def __init__(self,T_hot,T_cold,W_net=-1.,h_fg=100.):
        self.T_hot = T_hot
        self.T_cold = T_cold
        self.h_fg = h_fg
        self.W_net = W_net
        self.cop = 1./(T_hot / T_cold - 1.)
        self.Q_cold = self.W_net * self.cop
        self.Q_hot = -self.Q_cold - self.W_net
        self.m_dot = self.Q_hot / h_fg
        
    def __repr__(self):
        names = """T_hot T_cold m_dot h_fg Q_hot cop W_net Q_cold""".split()
        result = []
        for i in names:
            result.append(self.__getattribute__(i))
        return tabulate.tabulate(zip(names,result))
    
    def getCondenserStream(self):
        return HRHX_integral_model.streamExample2(
            self.h_fg,self.T_hot,-self.m_dot,1.0,self.h_fg)
    def getEvaporatorStream(self):
        deltaH = - self.Q_cold / self.Q_hot * self.h_fg
        return HRHX_integral_model.streamExample2(
            0, self.T_cold, -self.m_dot, 1.0, deltaH)

class CarnotPair(object):
    def __init__(self,T_hot,T_reject_power,T_reject_cooler,T_cold,m_dot=1):
        self.ce = CarnotEngine(T_hot,T_reject_power,m_dot)
        self.cc = CarnotCooler(T_reject_cooler, T_cold, -self.ce.W_net)
    def __repr__(self):
        return """Engine
{}
Cooler
{}""".format(self.ce.__repr__(), self.cc.__repr__())
    def display(self):
        import matplotlib.pyplot as plt
        plt.figure()
        plt.title("Internal streams")
        streams = [self.ce.getBoilerStream(), self.ce.getCondenserStream(),
                   self.cc.getCondenserStream(), self.cc.getEvaporatorStream()]
        QQ = [self.ce.Q_hot, self.ce.Q_cold, self.cc.Q_hot, self.cc.Q_cold]
        for stream, Q in zip(streams,QQ):
            plotStream(stream,(-0.1 * Q,1.1*Q))

class CarnotSystem(object):
    names = ["Boiler","Power condenser","Cooler condenser","Evaporator"]
        
    def __init__(self, carnot_pair, hot_stream, reject_power_stream,
                 reject_cooler_stream, cold_stream):
        self.carnot_pair = carnot_pair
        self.hot_stream = hot_stream
        self.reject_power_stream = reject_power_stream
        self.reject_cooler_stream = reject_cooler_stream
        self.cold_stream = cold_stream
        
        self.hotHX = HRHX_integral_model.counterflow_integrator(
            self.carnot_pair.ce.getBoilerStream(),self.hot_stream)
        self.engRejHX = HRHX_integral_model.counterflow_integrator(
            self.reject_power_stream,self.carnot_pair.ce.getCondenserStream())
        self.coolRejHX = HRHX_integral_model.counterflow_integrator(
            self.reject_cooler_stream,self.carnot_pair.cc.getCondenserStream())
        self.coldHX = HRHX_integral_model.counterflow_integrator(
            self.carnot_pair.cc.getEvaporatorStream(), self.cold_stream)
            
        self.HXs = [self.hotHX, self.engRejHX, self.coolRejHX, self.coldHX]
        self.QQ = [self.carnot_pair.ce.Q_hot, -self.carnot_pair.ce.Q_cold,
              -self.carnot_pair.cc.Q_hot, self.carnot_pair.cc.Q_cold]
                    
    def calcUA(self,output=False):
        UAs = []
        for HX, Q, name in zip(self.HXs, self.QQ, self.names):
            try:
                UA,eff = HX.calcUA(Q,True)
                UAs.append(UA)
                if output:
                    plt.figure()
                    HRHX_integral_model.plotFlow(HX,None,Q)
                    plt.title(name)
                    plt.xlabel("Heat flow (kW)")
                    plt.ylabel("Temperature (K)")
            except StandardError as e:
                print name
                print e
        if output:
            plt.figure()
            plt.bar(range(4),UAs,tick_label=self.names)
            plt.ylabel("UA kW/K")
            #plt.gca().set_yscale('log')
            print tabulate.tabulate(zip(self.names,UAs),["Component","UA (kW/K)"])
            print
        return np.array(UAs)
        
    def calcDeltaT(self,output=False):
        DeltaT = []
        for HX, Q, name in zip(self.HXs, self.QQ, self.names):
            dT = HX.calcDistanceT(Q)
            DeltaT.append(dT)
        if output:
            plt.figure()
            plt.bar(range(len(self.names)),DeltaT,tick_label=self.names)
            plt.ylabel("min DeltaT")
            #plt.gca().set_yscale('log')
            print tabulate.tabulate(zip(self.names,DeltaT),headers=["Component","DeltaT (K)"])
            print
        return np.array(DeltaT)

    def calcStuff(self,text=False,plot=False):
        DeltaT = []
        Eff = []
        UA = []
        for HX, Q, name in zip(self.HXs, self.QQ, self.names):
            dt,eff,ua = HX.calcUA2(Q)
            DeltaT.append(dt)
            Eff.append(eff)
            UA.append(ua)
            
            if plot:
                plt.figure()
                HRHX_integral_model.plotFlow(HX,None,Q)
                plt.title(name)
        if plot:
            n = len(self.names)
            plt.figure()
            plt.subplot(3,1,1)
            plt.ylabel("DeltaT")
            plt.bar(range(n),DeltaT,tick_label=self.names)

            plt.subplot(3,1,2)
            plt.ylabel("Eff.")
            plt.bar(range(n),Eff)
            
            plt.subplot(3,1,3)
            plt.ylabel("UA")
            plt.bar(range(n),UA)
            
        if text:
            print tabulate.tabulate(zip(self.names,DeltaT,Eff,UA),
                                    headers=["Component","DeltaT (K)","Eff","UA"])
            print
        return np.array(DeltaT),np.array(Eff),np.array(UA)
        
def plotStream(stream, qrange):
    # First let's plot T(q)
    q1 = np.linspace(*qrange)
    T1 = stream.T(q1)
    l=plt.plot(q1,T1,'-')[0]

def main(x=(400.,320.,320.,275.,1.),text=False,plot=False,collapse=True):
    if text:
        print "x = ", x
    cp = CarnotPair(*x)
    #print cp
    #cp.display()
    cs = CarnotSystem(carnot_pair=cp,
                      hot_stream=HRHX_integral_model.streamExample1(500,2),
                      reject_power_stream=HRHX_integral_model.streamExample1(313,15),
                      reject_cooler_stream=HRHX_integral_model.streamExample1(313,25),
                      cold_stream=HRHX_integral_model.streamExample1(279,50))
    objective = cs.carnot_pair.cc.Q_cold
    DeltaT,Eff,UA = cs.calcStuff(text=text,plot=plot)
    
    if collapse:
        return np.concatenate([DeltaT,[objective],UA])
    else:
        return DeltaT,objective,UA

def jacobian(f,x0,args=[],epsilon=0.001):
    oneOnEps = 1./epsilon
    n = len(x0)
    f0 = f(x0,*args)
    shape = (len(f0),len(x0))
    result = np.zeros(shape)
    for i in range(n):
        dx = np.zeros(n)
        dx[i] = epsilon
        x = x0 + dx
        fx = f(x,*args)
        df = fx - f0
        dfdx = df * oneOnEps
        result[:,i] = dfdx
    return result

class problem(object):
    def __init__(self,x0):
        self.calls=[0]
        self.x = [x0]
        self.DeltaT,self.y,self.UA = main(x0,collapse=False)
        self.UA0 = self.cost()
        self.constraints = []
        for i in range(6):
            self.constraints.append(dict(type="ineq",
                                         fun=lambda(x):self.evalConstraints(x)[i]))
        self.i=1
    def objective(self,x):
        print self.i," Objective at x = ", x
        self.i+=1
        self.calls.append(-1)
        self.x.append(x)
        self.DeltaT,self.y,self.UA = main(x,collapse=False)
        if any(self.DeltaT <= 0):
            return 0
        else:
            return -self.y
    def cost(self):
        return sum(self.UA)
    def evalConstraints(self,x):
        #print self.i, " Constraints at x = ", x
        self.i+=1
        self.calls.append(1)
        self.x.append(x)
        self.DeltaT,self.y,self.UA = main(x,collapse=False)
        
        result = []
        result.append(x[4]-0.1)
        for i in range(4):
            result.append(self.DeltaT[i]-0.1)
        result.append(self.UA0 - self.cost())
        
        return result
        
    def minimize(self):
        import scipy.optimize
        opt = scipy.optimize.minimize(self.objective,
                                    self.x[-1],
                                    constraints=self.constraints,
                                    options={"maxiter":100,"disp":True},
                                    method="COBYLA")
                                    #bounds=[(250,None)]*4+[(0.1,2)])
        return opt

if __name__ == "__main__":
    x0=np.array((400.,320.,320.,275.,1.))
    
    DeltaT0,y0,UA0=main(x0,text=True,plot=True,collapse=False)
    print "DeltaT =", DeltaT0
    print "objective = {} kW".format(y0)
    print "UAs = {}, sum = {}".format(UA0,sum(UA0))
    
    if False:
        j = jacobian(main,x0,[False,False])
        print j
        print "Rank of jacobian:",np.linalg.matrix_rank(j)
        plt.figure()
        plt.subplot(3,2,1)
        plt.gca().matshow(j[0:4,0:4])
        plt.subplot(3,2,2)
        plt.gca().matshow(j[0:4,4:5])
        plt.subplot(3,2,3)
        plt.gca().matshow(j[4:5,0:4])
        plt.subplot(3,2,4)
        plt.gca().matshow(j[4:5,4:5])
        plt.subplot(3,2,5)
        plt.gca().matshow(j[5:9,0:4])
        plt.subplot(3,2,6)
        plt.gca().matshow(j[5:9,4:5])
        
        # Solution
        # 1. Satisfy the constraint by travelling along the given path
        # 2. Optimize the objective by travelling in the null space of the objective.
        costm = np.array([1,1,1,1,0])
    print
    p=problem(x0)
    print
    opt=p.minimize()
    x1 = opt.x
    
    DeltaT1,y1,UA1=main(x1,text=True,plot=True,collapse=False)
    print "DeltaT =", DeltaT1
    print "y = {} kW".format(y1)
    print "UA = {}, sum = {}".format(UA1,sum(UA1))
    
    plt.figure()
    plt.bar(np.arange(4),UA0,0.3,tick_label=CarnotSystem.names)
    plt.bar(np.arange(4)+0.3,UA1,0.3,color="red")
        