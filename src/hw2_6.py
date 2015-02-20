# -*- coding: utf-8 -*-
"""
Created on Mon Feb 09 15:00:58 2015

@author: nfette
"""

from CoolProp.CoolProp import PropsSI
import numpy as np
from hw2_1 import *
import matplotlib.pyplot as plt
from scipy.special import logit, expit
from scipy.optimize import minimize

import inspect
me = inspect.getfile(inspect.currentframe())

def zeroOrderModel(T0, T1, T2, R0, R1, R2, args):
        # Initial conditions
        T0i, T2i = map(CelsiusToKelvin, (T0 - 10, T2 - 10))
        # Solver loop
        for i in range(args):
            T1i = 0.5 * (T0i+T2i)
            # Calculate the heat flows due to heat transfer
            dt0, dt1, dt2 = (T0-T0i, T1-T1i, T2-T2i)
            Q0, Q1, Q2 = dt0 / R0, dt1 / R1, dt2 / R2
            # Calculate the approximate COP
            COP = COP_cooling_reversible(T0i, T1i, T2i)
    
            f,g,h = COP_cooling_partial_Tei(T0i, T1i, T2i), \
                COP_cooling_partial_Tci(T0i, T1i, T2i), \
                COP_cooling_partial_Thi(T0i, T1i, T2i)
            a = f + 0.5 * g
            b = h + 0.5 * g
            
            A = np.zeros((2,2))
            A[0,0] = (-1. / R0) - (a/R2)*dt2
            A[0,1] = (-b/R2) * dt2 + (COP / R2)
            A[1,0] = -0.5 / R1 + (a/R2)*dt2
            A[1,1] = -0.5 / R1 + (b/R2)*dt2 - (1/R2)*(1+COP)
            B = np.ones(2)
            B[0] = (-1 / R0) * dt0 + (COP / R2) * dt2
            B[1] = (-1 / R1) * dt1 - (1 + COP) / R2 * dt2
            DT = np.linalg.solve(A,B)
            # update
            T0i = T0i + DT[0]
            T2i = T2i + DT[1]
                        
        return COP, Q0, Q1, Q2, T0i, T1i, T2i

T0, T1, T2 = map(CelsiusToKelvin, (-20, 50, 200))
U = 500.
Atot = 10.

def objective1(A012, args):
    A0, A1, A2 = A012
    R0, R1, R2 = 1 / (U * A0), 1 / (U * A1), 1 / (U * A2)

    COP, Q0, Q1, Q2, T0i, T1i, T2i = zeroOrderModel(T0, T1, T2, R0, R1, R2, args)
    return -Q0

cons = ({'type': 'ineq',
         'fun' : lambda x: np.array([x[0]]),
         'jac' : lambda x: np.array([1.,0.,0.])},
        {'type': 'ineq',
         'fun' : lambda x: np.array([Atot - x[0]]),
         'jac' : lambda x: np.array([-1.0, 0., 0.])},
        {'type': 'ineq',
         'fun' : lambda x: np.array([x[1]]),
         'jac' : lambda x: np.array([0.0, 1., 0.])},
        {'type': 'ineq',
         'fun' : lambda x: np.array([Atot - x[1]]),
         'jac' : lambda x: np.array([0., -1., 0.])},
        {'type': 'ineq',
         'fun' : lambda x: np.array([x[2]]),
         'jac' : lambda x: np.array([0., 0., 1.])},
        {'type': 'ineq',
         'fun' : lambda x: np.array([Atot - x[2]]),
         'jac' : lambda x: np.array([0., 0., -1.])},
        {'type': 'eq',
         'fun' : lambda x: np.array([x[0] + x[1] + x[2] - Atot]),
         'jac' : lambda x: np.array([1., 1., 1.])})

def min1():
    iters = 4
    return minimize(objective1, (2.,3.,5.), constraints = cons,
                       method = 'SLSQP', args=(iters,))

def min2(result):
    iters = 6
    return minimize(objective1, result.x, constraints = cons,
                       method = 'SLSQP', args=(iters,))
    
def main():
    # Obtain the baseline configuration from example 2.2
    result = min1()
    print result
    print
    result2 = min2(result)
    print result2
    A0, A1, A2 = result2.x
    UA0, UA1, UA2 = U * A0, U * A1, U * A2
    R0, R1, R2 = 1 / UA0, 1 / UA1, 1 / UA2
    COP, Q0, Q1, Q2, T0i, T1i, T2i = zeroOrderModel(T0, T1, T2, R0, R1, R2, 6)
    UA0_nominal, COP_nominal, Q0_nominal = U * A0, COP, Q0
    T0i_nom, T1i_nom, T2i_nom = T0i, T1i, T2i
    print
    print "A0, A1, A2 = {}, {}, {} m^2".format(A0, A1, A2)
    print "T = {},{},{},{},{},{} K".format(T0i, T0, T1, T1i, T2i, T2)
    print "Q = {}, {}, {} W".format(Q0, Q1, Q2)
    print "COP = {}".format(COP)
    print "Residuals {:e}, {:e} W".format(Q0 + Q1 + Q2, COP * Q2 - Q0)
    
    plt.close('all')
    
    fig,ax1 = plt.subplots()
    ax2 = ax1.twinx()
    x = np.array([0,1,2])
    width = 0.33
    ax1.bar(x, [Q0, -Q1, Q2], width, color='red')
    ax2.bar(x+width, [UA0, UA1, UA2], width, color='blue')
    ax1.set_ylabel('Unsigned Heat flows',color='red')
    ax1.set_xticks(x + width)
    ax1.set_xticklabels(('Cold','Mid','Hot'))
    ax2.set_ylabel('Exchanger UA',color='blue')
    plt.savefig('{}.figure{}.png'.format(me,plt.gcf().number))
    
    # Vary the evaporator heat exchanger.
    A_evap = np.linspace(1,10)
    a_UA = []
    a_Q0, a_Q2 = [], []
    a_COP = []
    a_T0i, a_T1i, a_T2i = [], [], []
    for A0_new in A_evap:
        UA_new = U * A0_new
        a_UA.append(UA_new)
        COP, Q0, Q1, Q2, T0i, T1i, T2i = \
            zeroOrderModel(T0, T1, T2, 1/UA_new, R1, R2, 6)
        a_Q0.append(Q0)
        a_Q2.append(Q2)
        a_COP.append(COP)
        a_T0i.append(T0i)
        a_T1i.append(T1i)
        a_T2i.append(T2i)
    
    # Vary the middle temperature heat exchanger.
    A_cond = np.linspace(1,10)
    a2_UA = []
    a2_Q0, a2_Q2 = [], []
    a2_COP = []
    a2_T0i, a2_T1i, a2_T2i = [], [], []
    for A1 in A_cond:
        UA_new = U * A1
        a2_UA.append(UA_new)
        COP, Q0, Q1, Q2, T0i, T1i, T2i = zeroOrderModel(T0, T1, T2, R0, 1/UA_new, R2, 6)
        a2_Q0.append(Q0)
        a2_Q2.append(Q2)
        a2_COP.append(COP)
        a2_T0i.append(T0i)
        a2_T1i.append(T1i)
        a2_T2i.append(T2i)
    
    fig1, ax1 = plt.subplots()
    ax1.plot(UA0_nominal, Q0_nominal, 'o')
    ax2 = ax1.twinx()
    ax2.plot(UA0_nominal, COP_nominal, 'ro')
    ax1.plot(a_UA, a_Q0, 'b')

    #ax1.plot(a_UA, a_Q2, 'g')
    ax1.set_xlabel('Evaporator UA (W/K)')
    ax1.set_ylabel('Capacity Q0 (W)', color='b')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')
    
    ax2.plot(a_UA, a_COP, 'r')
    ax2.set_ylabel('COP = Q0 / Q2', color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')

    plt.savefig('{}.figure{}.png'.format(me,plt.gcf().number))
    
    plt.figure()
    plt.plot(a_UA, a_T0i, a_UA, a_T1i, a_UA, a_T2i)
    plt.gca().axhline(T0,ls='--',color='b')
    plt.gca().axhline(T1,ls='--',color='g')
    plt.gca().axhline(T2,ls='--',color='r')
    plt.gca().axvline(UA0_nominal,color='k')
    plt.xlabel('Evaporator UA (W/K)')
    plt.ylabel('Temperature (K)')
    plt.savefig('{}.figure{}.png'.format(me,plt.gcf().number))
    # plt.show()

    return

if __name__ == "__main__":
    main()
