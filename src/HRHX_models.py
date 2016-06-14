# -*- coding: utf-8 -*-
"""
HRHX_models.py
2016-02-29 by Nicholas Fette and Andrew Hickey
Models for the heat recovery heat exchanger.
"""
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve
from scipy.integrate import quad
import scipy.optimize as op

class HRHX_model1():
    """An effectiveness model. Assumes and requires that exhaust stream capacity
    (C = mdot * C_p) is less than that of the HTF stream.
    Uses an estimate of mean heat capacity, treated as constant.
    Ignores the possibility of interior pinch.
    
    Usage
    -----
    e,T,P,m = 0.96, 270., 100, 0.80
    HRHX = HRHX_model(e,T,P,m)
    T_in, m_in, htf = 50, 2, HTFLookup(filename)
    Q, T_HTF_outlet, T_exhaust_outlet, DeltaP_HTF, DeltaP_exhaust \\
        = HRHX(T_in, m_in, htf)
    """
    def __init__(self, effectiveness, T_exhaust_inlet, P_exhaust_inlet, m_exhaust):
        """Get things going."""
        self.effectiveness = effectiveness # [fraction], 0 to 1
        self.T_exhaust_inlet = T_exhaust_inlet # [C]
        self.P_exhaust_inlet = P_exhaust_inlet # [Pa]
        self.m_exhaust = m_exhaust # [kg/s]
    def __call__(self, T_HTF_inlet, m_in, htf):
        """Compute the outputs for these inputs."""
        # Function maps inputs ( T [K], P [Pa] ) to outputs (C [J/kg-K])
        # We apply it:  inputs ( T [C], P [Pa] ) to outputs (C [J/kg-K])
        T_min, T_max = T_HTF_inlet, self.T_exhaust_inlet
        T_mean = 0.5 * (T_min + T_max)
        C_p_exhaust = PropsSI('C',
                                     'T', T_mean + 273.15,
                                     'P', self.P_exhaust_inlet,
                                     'Air')
        C_exhaust = C_p_exhaust * self.m_exhaust  # [J/K]
        C_HTF = htf.lookup('C',T_mean) * m_in

        # This should be true
        C_min, C_max = C_exhaust, C_HTF
        C_r = C_min / C_max
        if C_r > 1:
            raise ValueError("HTF stream should have greater capacity (C_r = {})".format(C_r))
        DeltaT_max = self.T_exhaust_inlet - T_HTF_inlet
        Q_max = C_min * DeltaT_max
        Q = self.effectiveness * Q_max
        T_HTF_outlet = T_HTF_inlet + Q / C_HTF
        T_exhaust_outlet = self.T_exhaust_inlet - Q / C_exhaust

        # TODO: adjust for (water) flow rate
        DeltaP_HTF = -6e3 # [Pa]
        DeltaP_exhaust = -300 # [Pa]
        
        return Q, T_HTF_outlet, T_exhaust_outlet, DeltaP_HTF, DeltaP_exhaust

HRHX1_default = HRHX_model1(0.96, 270., 100, 0.80)

class HRHX_model2():
    """An effectiveness model. Assumes and requires that exhaust stream capacity
    (C = mdot * C_p) is less than that of the HTF stream.
    Uses enthalpy for calculations instead of average C_p, but ignore
    the possibility of interior pinch.
    
    Usage
    -----
    e,T,P,m = 0.96, 270., 100, 0.80
    HRHX = HRHX_model(e,T,P,m)
    T_in, m_in, htf = 50, 2, HTFLookup(filename)
    Q, T_HTF_outlet, T_exhaust_outlet, DeltaP_HTF, DeltaP_exhaust \\
      = HRHX(T_in, m_in, htf)
        
    Args:
        effectiveness (float): a ratio between 0 and 1
        T_exhaust_inlet (float): temperature in [C]
        P_exhaust_inlet (float): pressure in [Pa]
        m_exhaust (float): a mass flow rate in [kg/s]
    """
    def __init__(self, effectiveness, T_exhaust_inlet, P_exhaust_inlet, m_exhaust):
        
        self.effectiveness = effectiveness # [fraction], 0 to 1
        self.T_exhaust_inlet = T_exhaust_inlet # [C]
        self.P_exhaust_inlet = P_exhaust_inlet # [Pa]
        self.m_exhaust = m_exhaust # [kg/s]
    def __call__(self, T_HTF_inlet, m_in, htf):
        """Apply the model to a given HTF inlet condition.
        
        Args:
            T_HTF_inlet(float)
            m_in(float)
            htf(HTFLookup)
        
        Returns:
            Q
            T_HTF_outlet
            T_exhaust_outlet
            DeltaP_HTF
            DeltaP_exhaust
            
        """
        T_min, T_max = T_HTF_inlet, self.T_exhaust_inlet
        # Function maps inputs ( T [K], P [Pa] ) to outputs (H [J/kg])
        # We apply it:  inputs ( T [C], P [Pa] ) to outputs (H [J/kg])
        H_exhaust = lambda(T): PropsSI('H',
                                     'T', T + 273.15,
                                     'P', self.P_exhaust_inlet,
                                     'Air')
        H_exhaust_max = H_exhaust(T_max)
        DeltaH_exhaust_max = H_exhaust_max - H_exhaust(T_min)
        # Q_max = C_min * DeltaT_max
        Q_max = self.m_exhaust * DeltaH_exhaust_max
        Q = self.effectiveness * Q_max
        # Now find the temperatures that give this value.
        # For the exhaust stream we have an inverse.
        H_exhaust_outlet = H_exhaust_max  - DeltaH_exhaust_max * self.effectiveness
        # Function maps inputs ( T [K], P [Pa] ) to outputs (H [J/kg])
        # We apply it:  inputs ( T [C], P [Pa] ) to outputs (H [J/kg])
        T_exhaust_outlet = -273.15 + PropsSI('T',
                                            'H', H_exhaust_outlet,
                                            'P', self.P_exhaust_inlet,
                                            'Air')
        DeltaH_HTF = Q / m_in
        H_HTF_min = htf.h(T_HTF_inlet)
        self.H_HTF_max = H_HTF_min + DeltaH_HTF
        self.H_HTF_err = lambda(T): htf.h(T) - self.H_HTF_max
        self.T_HTF_outlet = fsolve(self.H_HTF_err, T_HTF_inlet)
        T_HTF_outlet = self.T_HTF_outlet[0]

        # TODO: adjust for (water) flow rate
        DeltaP_HTF = -6e3 # [Pa]
        DeltaP_exhaust = -300 # [kPa]
        
        return Q, T_HTF_outlet, T_exhaust_outlet, DeltaP_HTF, DeltaP_exhaust

class HRHX_model3(object):
    """An effectiveness model. Assumes and requires that exhaust stream capacity
    (C = mdot * C_p) is less than that of the HTF stream.
    Uses enthalpy for calculations instead of average C_p, but ignore
    the possibility of interior pinch.
    
    Usage
    -----
    e,T,P,m = 0.96, 270., 100, 0.80
    HRHX = HRHX_model(e,T,P,m)
    T_in, m_in, htf = 50, 2, HTFLookup(filename)
    Q, T_HTF_outlet, T_exhaust_outlet, DeltaP_HTF, DeltaP_exhaust \\
      = HRHX(T_in, m_in, htf)
        
    Args:
        effectiveness (float): a ratio between 0 and 1
        T_exhaust_inlet (float): temperature in [C]
        P_exhaust_inlet (float): pressure in [Pa]
        m_exhaust (float): a mass flow rate in [kg/s]
    """
    def __init__(self,UA):
        self.UA=UA
    def __call__(self,T_cold,T_hot):
        Q=0
        out1=0
        out2=0
        return Q,out1,out2

def calcUA(T_cold, T_hot, Q):
    """Compute UA value given the heat flux.
    Based on equation 10 in prospectus.
    
    Args
    ----
        T_cold (function): maps cumulative heat flux relative to cold inlet
                           to temperature of cold stream (°C).
        T_hot (function): maps cumulative heat flux relative to hot inlet
                          to temperature of hot stream (°C).
        Q (float): the total heat flux (W).
    
    Returns
    -------
        UA (float): the overall heat exchanger coefficient (W/°C)
    """
    integrand = lambda(q): 1 / (T_hot(Q-q)-T_cold(q))
    result = quad(integrand,0,Q)
    return result

def calcQ(T_cold,T_hot,UA,Q_max):
    """Compute heat flux given the UA value.
    Based on equation 10 in prospectus. Searches over Q from 0 to Q_max.
    
    Args
    ----
        T_cold (function): maps cumulative heat flux relative to cold inlet
                           to temperature of cold stream (°C).
        T_hot (function): maps cumulative heat flux relative to hot inlet
                          to temperature of hot stream (°C).
        UA (float): the overall heat exchanger coefficient (W/°C)
        Q_max (float): the maximum heat flux to search over (W).
    
    Returns
    -------
        Q (float): the total heat flux (W).
    """
    f = lambda(Q):calcUA(T_cold,T_hot,Q)[0]-UA
    result = fsolve(f,0)
    return result
    
def calcQmax(q_cold,q_hot,T_min,T_max):
    T_guess = 0.5 * (T_min + T_max)
    f = lambda(t):q_cold(t)+q_hot(t)

    cons = ({'type': 'ineq', 'fun': lambda x:  x-T_min},     # >= 0
            {'type': 'ineq', 'fun': lambda x: T_max-x})
    result = op.minimize(f,T_guess,constraints=cons,options={'disp':True})
    #result = bracket(f,T_guess,bounds=(T_min,T_max))
    return result.fun[0]

if __name__ == "__main__":
    T_cold=lambda(q):q
    T_hot=lambda(q):1-q
    import numpy as np
    import matplotlib.pyplot as plt
    q=np.linspace(0,1)
    Q = 0
    tc = T_cold(q)
    th = T_hot(q)
    #plt.plot(q,tc,Q-q,th)
    Q = 0.9
    plt.figure(1)
    plt.plot(q,tc,Q-q,th)
    UA_expect = Q / (1-Q)
    UA = calcUA(T_cold,T_hot,Q)[0]
    print UA_expect, UA
    
    calcQ(T_cold,T_hot,9,None)
    QQ = 1-np.logspace(-3,0)
    UA = np.zeros_like(QQ)
    for i in range(len(QQ)):
        UA[i] = calcUA(T_cold, T_hot, QQ[i])[0]
    plt.figure(2)
    plt.loglog(1-QQ,UA)
    
    T_cold = lambda(q):q
    T_hot = lambda(q):1-q**2
    q_cold = lambda(t):t
    q_hot = lambda(t):np.sqrt(1-t)
    
#    T_cold = lambda(q):q
#    T_hot = lambda(q):(q-1)**2
#    q_cold = lambda(t):t
#    q_hot = lambda(t):1-np.sqrt(t)
    
    T_cold = lambda(q):0*q
    T_hot = lambda(q):1-q**2
    q_cold = np.vectorize(lambda(t):0 if t <= 0 else 2)
    q_hot = lambda(t):np.sqrt(1.-t)
    
    tc = T_cold(q)
    th = T_hot(q)
    plt.figure(3)
    Q_max = calcQmax(q_cold,q_hot,0,1)
    Q=Q_max
    plt.plot(q,tc,Q-q,th)
    
    QQ = Q_max*(1-np.logspace(-5,0))
    UA = np.zeros_like(QQ)
    for i in range(len(QQ)):
        UA[i] = calcUA(T_cold, T_hot, QQ[i])[0]
    plt.figure(2); plt.cla()
    plt.loglog(Q_max-QQ,UA)
    