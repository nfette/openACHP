# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 15:48:09 2016

@author: nfette
"""
import tabulate

class CarnotEngine(object):
    """Temperature in K"""
    def __init__(self,T_hot,T_cold,m_dot=1.,h_fg=1.):
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

class CarnotCooler(object):
    def __init__(self,T_hot,T_cold,W_net=-1.,h_fg=1.):
        self.T_hot = T_hot
        self.T_cold = T_cold
        self.h_fg = h_fg
        self.W_net = W_net
        self.cop = 1./(T_hot / T_cold - 1.)
        self.Q_cold = - self.W_net * self.cop
        self.Q_hot = -self.Q_cold - self.W_net
        self.m_dot = self.Q_hot / h_fg
        
    def __repr__(self):
        names = """T_hot T_cold m_dot h_fg Q_hot cop W_net Q_cold""".split()
        result = []
        for i in names:
            result.append(self.__getattribute__(i))
        return tabulate.tabulate(zip(names,result))
    
    def getCondenserStream(self):
        pass
    def getEvaporatorStream(self):
        pass

class CarnotPair(object):
    def __init__(self,T_hot,T_reject_power,T_reject_cooler,T_cold,m_dot=1):
        self.ce = CarnotEngine(T_hot,T_reject_power)
        self.cc = CarnotCooler(T_reject_cooler, T_cold, -ce.W_net)
    def __repr__(self):
        return """Engine
{}
Cooler
{}""".format(self.ce.__repr__(), self.cc.__repr__())
        
cp = CarnotPair(400.,313.,313.,279.)
print cp