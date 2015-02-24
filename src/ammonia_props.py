# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 11:36:55 2015

@author: nfette
"""

from __future__ import print_function
import ees_interface
from collections import namedtuple
import sys
thismodule = sys.modules[__name__]

molecular_mass_water = 18.015 # g/mol
molecular_mass_ammonia = 17.031 # g/mol

defaultPath = r"C:\EES32\Userlib\EES_System\nh3h2o.dlp"

availableCodes = [123, 128, 137, 138, 148, 158, 168, 178,\
    234, 235, 238, 248, 258, 268, 278]
availableCodeStringsForward = {}
availableCodeStringsBackward = {}

def splitCode1(code):
    return [code / 100, code / 10 % 10, code % 10]

def splitCode0(code):
    return [code / 100 - 1, code / 10 % 10 - 1, code % 10 - 1]

# Eg, ['T','P','x'],outvars -> 123
def reverseCodeLookup(propKeys,outvars):
    a=[outvars.index(i)+1 for i in propKeys]
    a.sort()
    code = int('{0[0]}{0[1]}{0[2]}'.format(a))
    keysSorted = [outvars[i-1] for i in a]
    return code, keysSorted

def standardOrder(state,code):
    return tuple([state[i] for i in availableCodeStringsForward[code]])

standardOutvars = ['T', 'P', 'x', 'h', 's', 'u', 'v', 'Qu']
State = namedtuple('State',standardOutvars)
def gibbs(self):
    return self.h - self.T * self.s
State.gibbs = gibbs
def molarMass(self):
    return 1. / (self.x / molecular_mass_ammonia \
        + (1. - self.x) / molecular_mass_water)
State.molarMass = molarMass
dgdxvars = ['dgdx','mu1','mu2']
State2 = namedtuple('State2',dgdxvars)
State3 = namedtuple('State3',['dhdx','h1','h2'])

class ammoniaWaterFunc:
    def __init__(self,mydll,code,S=""):
        self.mydll = mydll
        if type(code) == int:
            self.code = code
        elif type(code) == str:
            self.code = availableCodeStringsBackward[code]
        self.S = S
    def getCallFormat(self,*args):
        callFormat, _, outvars = self.mydll.getCallFormat(self.S,(self.code,)+args)
        invars = availableCodeStringsForward[self.code]
        return callFormat, invars, outvars
    def getInputUnits(self,inarglist=[0]):
        return self.mydll.getInputUnits(self.S,[self.code]+inarglist)[1:]
    def getOutputUnits(self,inarglist=[0]):
        return self.mydll.getOutputUnits(self.S,[self.code]+inarglist)
    def call(self,*args):
        return State(*self.mydll.call(self.S,(self.code,)+args)[1])
    def dgdxetc(self,deltax=0.0001,**args):
        state0 = self.call(*standardOrder(args,self.code))
        args['x'] += deltax
        state1 = self.call(*standardOrder(args,self.code))
        dgdx = (state1.gibbs() - state0.gibbs()) / deltax
        mu1 = state0.gibbs() + (1 - state0.x) * dgdx
        mu2 = state0.gibbs() - state0.x * dgdx
        return State2(dgdx,mu1,mu2)
    def dhdxetc(self,deltax=0.0001,**args):
        state0 = self.call(*standardOrder(args,self.code))
        args['x'] += deltax
        state1 = self.call(*standardOrder(args,self.code))
        dhdx = (state1.h - state0.h) / deltax
        h1 = state0.h + (1 - state0.x) * dhdx
        h2 = state0.h - state0.x * dhdx
        return State3(dhdx,h1,h2)

class AmmoniaProps:
    def __init__(self,path = defaultPath):
        self.mydll = ees_interface.EES_DLP(path)
        _,_,self.outvars = self.mydll.getCallFormat()
        #thismodule.State = namedtuple('State',self.outvars)
        #State.gibbs = gibbs
        for code in availableCodes:
            inputs = [self.outvars[i] for i in splitCode0(code)]
            availableCodeStringsForward[code] = inputs
            availableCodeStringsBackward[''.join(inputs)] = code
    def props(self,code):
        return ammoniaWaterFunc(self.mydll,code)

if __name__ == "__main__":
    myprops = AmmoniaProps(defaultPath)
    f1 = myprops.props(123)
    callFormat, invars, outvars = f1.getCallFormat()
    print(callFormat)
    print("Invars, outvars:{},{}".format(invars,outvars))
    print(f1.getInputUnits())
    print(f1.getOutputUnits())

    print(myprops.props(123).call(450, 10, 0.5))
    print(myprops.props('TPx').call(450, 10, 0.5))
    try:
        print(myprops.props('TPu').call(450, 10, 0.5))
    except Exception as e:
        print(e.__repr__())
        
