# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 11:36:55 2015

@author: nfette
"""

import ees_interface
from collections import namedtuple
import sys
import numpy as np
import tabulate
import json
thismodule = sys.modules[__name__]

molecular_mass_water = 18.015 # g/mol
molecular_mass_ammonia = 17.031 # g/mol

defaultPath = r"C:\EES32\Userlib\EES_System\nh3h2o.dlp"

availableCodes = [123, 128, 137, 138, 148, 158, 168, 178,\
    234, 235, 238, 248, 258, 268, 278]
availableCodeStringsForward = {}
availableCodeStringsBackward = {}

standardOutvars = ['T', 'P', 'x', 'h', 's', 'u', 'v', 'Qu']
standardUnits = ['K', 'bar', ' ', 'kJ/kg', 'kJ/kg-K', 'kJ/kg', 'm^3/kg', ' ']
State = namedtuple('State',standardOutvars)
StateType = np.dtype(dict(names=standardOutvars,formats='f'*8))

def convert_state_list_to_array(state_list):
    """Input a list of State objects, ouput an array of StateType."""
    state_array = np.zeros(len(state_list), dtype=StateType)
    for i, s in enumerate(state_list):
        state_array[i] = s
    return state_array

class CStateTable:
    def __init__(self, states, labels=None):
        """Input states, an array of StateType, and labels, a list of
        labels for the points. Makes printing table convenient."""
        self.states = states
        self.labels = labels
        if labels is None:
            self.labels = [''] * len(states)
    def doit(self, **kwargs):
        return tabulate.tabulate([[p] \
                                  + [s[name] for name in StateType.names]
                                  for p, s in zip(self.labels, self.states)],
                                 StateType.names,
                                 **kwargs)
    def __repr__(self):
        return self.doit()
    def _repr_html_(self):
        return self.doit(tablefmt="html")

    def toJSON(self):
        output = {}
        for p, s in zip(self.labels, self.states):
            node = {}
            for key in StateType.names:
                node[key] = s[key].astype(float)
            output[p] = node
        return json.dumps(output)

def splitCode1(code):
    return [code // 100, code // 10 % 10, code % 10]

def splitCode0(code):
    return (code // 100 - 1, code // 10 % 10 - 1, code % 10 - 1)
for code in availableCodes:
    inputs = [standardOutvars[i] for i in splitCode0(code)]
    availableCodeStringsForward[code] = inputs
    availableCodeStringsBackward[''.join(inputs)] = code

# Eg, ['T','P','x'],outvars -> 123
def reverseCodeLookup(propKeys,outvars=standardOutvars):
    a=[outvars.index(i)+1 for i in propKeys]
    a.sort()
    code = int('{0[0]}{0[1]}{0[2]}'.format(a))
    keysSorted = [outvars[i-1] for i in a]
    return code, keysSorted

# Eg, P=pval,s=sval,x=xval -> (235, pval, xval, sval)
def encode(**kwargs):
    code, keysSorted = reverseCodeLookup(kwargs)
    return (code,) + tuple([kwargs[key] for key in keysSorted])

def standardOrder(state,code):
    return tuple([state[i] for i in availableCodeStringsForward[code]])


def gibbs(self):
    return self.h - self.T * self.s
State.gibbs = gibbs
def molarMass(self):
    return 1. / (self.x / molecular_mass_ammonia \
        + (1. - self.x) / molecular_mass_water)
State.molarMass = molarMass
def isSubcooled(self):
    return self.Qu < 0
def isSuperheated(self):
    return self.Qu > 1
State.isSubcooled = isSubcooled
State.isSuperheated = isSuperheated

dgdxvars = ['dgdx','mu1','mu2']
State2 = namedtuple('State2',dgdxvars)
State3 = namedtuple('State3',['dhdx','h1','h2'])    

def massFractionToMolar(w):
    """Converts ammonia mass fraction w to ammonia molar fraction x."""
    return (w / molecular_mass_ammonia) \
        / (w / molecular_mass_ammonia + (1. - w) / molecular_mass_water);

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
    def __call__(self,*args):
        return State(*self.mydll.call(self.S,(self.code,)+args)[1])
    def dgdxetc(self,deltax=0.0001,**args):
        state0 = self(*standardOrder(args,self.code))
        args['x'] += deltax
        state1 = self(*standardOrder(args,self.code))
        dgdx = (state1.gibbs() - state0.gibbs()) / deltax
        mu1 = state0.gibbs() + (1 - state0.x) * dgdx
        mu2 = state0.gibbs() - state0.x * dgdx
        return State2(dgdx,mu1,mu2)
    def dhdxetc(self,deltax=0.0001,**args):
        state0 = self(*standardOrder(args,self.code))
        args['x'] += deltax
        state1 = self(*standardOrder(args,self.code))
        dhdx = (state1.h - state0.h) / deltax
        h1 = state0.h + (1 - state0.x) * dhdx
        h2 = state0.h - state0.x * dhdx
        return State3(dhdx,h1,h2)
    def massSpecificHeat(self,deltaT=1e-4,**args):
        print("deltaT = {} K".format(deltaT))
        state0 = self(*standardOrder(args,self.code))
        args['T'] += deltaT
        state1 = self(*standardOrder(args,self.code))
        dhdT = (state1.h - state0.h) / deltaT
        return dhdT
        
class AmmoniaProps:
    """
    Imports the EES NH3H2O library, if installed.
    """
    def __init__(self,path = defaultPath):
        """
        Args
        ----
            path : (string)
                The full path to the DLL file to load.
        """
        self.mydll = ees_interface.EES_DLP(path)
        _,_,self.outvars = self.mydll.getCallFormat()
        #thismodule.State = namedtuple('State',self.outvars)
        #State.gibbs = gibbs
        self.props2v = np.vectorize(self.props2)
    def props(self,code):
        """Returns an instance of ammoniaWaterFunc instance for the given set
        of input parameters. The resulting object can be called as a functor.
        
        Args
        ----
            code : (string or number)
                Defines the three input parameters.
                Can be one of the numbers in availableCodes, or a string
                concatenation of the corresponding names.
        """
        return ammoniaWaterFunc(self.mydll,code)
    def props2(self,**kwargs):
        """Returns the state corresponding to the given set of inputs and
        values. To determine the units, you will have to use props().
        
        kwargs
        ------
        Choose three in a combination matching the availableCodes.
            T : (float)
                Temperature (K)
            P : (float)
                Pressure (bar)
            x : (float)
                Ammonia mass fraction (kg/kg)
            h : (float)
                Mass specific enthalpy (kJ/kg)
            s : (float)
                Mass specific entropy (kJ/kg-K)
            u : (float)
                Mass specific internal energy (kJ/kg)
            v : (float)
                Mass specific volume (m^3/kg)
            Qu : (float)
                Vapor quality (kg/kg)
            out : (string)
                [Optional] Just return the named variable instead of a State.
        """
        try:
            out = kwargs['out']
            del kwargs['out']
        except:
            out = None
        args=encode(**kwargs)
        #if args[0] not in availableCodes:
        #    raise ValueError("Input code not implemented: {}".format(args[0]))
        s,vals = self.mydll.call("",args)
        if s:
            #print(vals)
            raise KeyError('DLL returned {}'.format(s))
        
        if out:
            return vals[standardOutvars.index(out)]
        else:
            return State(*vals)
    def T(self,**kwargs):
        return self.props2(**kwargs).T
    def P(self,**kwargs):
        return self.props2(**kwargs).P
    def x(self,**kwargs):
        return self.props2(**kwargs).x
    def h(self,**kwargs):
        return self.props2(**kwargs).h
    def s(self,**kwargs):
        return self.props2(**kwargs).s
    def u(self,**kwargs):
        return self.props2(**kwargs).u
    def Qu(self,**kwargs):
        return self.props2(**kwargs).Qu
    def equilibriumStates(self, P, z):
        """ Return the liquid state at P,z,Qu=0,
        and the corresponding vapor state at P,T,Qu=1
        """
        liquid = self.props2(P=P,x=z,Qu=0)
        vapor = self.props2(P=P,T=liquid.T,Qu=1)
        return liquid,vapor
    def equilibriumStates2(self, P, z_vapor):
        """ Return the vapor state at P,z_vapor,Qu=1,
        and the equilibrium liquid state at P,T,Qu=0
        """
        vapor = self.props2(P=P,x=z_vapor,Qu=1)
        liquid = self.props2(P=P,T=vapor.T,Qu=0)
        return liquid,vapor
    def equilibriumStates3(self, P, T):
        """ Return the liquid state at P,T,Qu=0,
        and the corresponding vapor state at P,T,Qu=1
        """
        liquid = self.props2(P=P,T=T,Qu=0)
        vapor = self.props2(P=P,T=T,Qu=1)
        return liquid,vapor

if __name__ == "__main__":
    myprops = AmmoniaProps(defaultPath)
    f1 = myprops.props(123)
    callFormat, invars, outvars = f1.getCallFormat()
    print(callFormat)
    print("Invars, outvars:{},{}".format(invars,outvars))
    print(f1.getInputUnits())
    print(f1.getOutputUnits())

    print(myprops.props(123)(450, 10, 0.5))
    print(myprops.props('TPx')(450, 10, 0.5))
    print(myprops.props2(T=450, P=10, x=0.5))
    try:
        print(myprops.props('TPu')(450, 10, 0.5))
    except KeyError as e:
        print(e.__repr__())
        
    try:
        print(myprops.props2(T=450, P=10, u=0.5))
    except KeyError as e:
        print(e.__repr__())
        print("Successfully caught error.")

