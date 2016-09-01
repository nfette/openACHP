# -*- coding: utf-8 -*-
"""
Created on Mon Aug 08 07:50:51 2016

@author: nfette
"""
import CoolProp
#import numpy as np
from numpy import linspace, exp, log
import matplotlib.pyplot as plt
from scipy.integrate import quad, odeint
from hw2_1 import CelsiusToKelvin as C2K, KelvinToCelsius as K2C

class AdsorptionChillerSpec(object):
    """Variables
    a_hex = area of adsorber heat exchanger (m^2)
    cpv = specific heat of the refrigerant vapor from the evaporator
    cw = specific heat for the refrigerant( or exhaust water)
    c1 = specific heat for the metallic adsorber (am)
    c2 = specific heat for the adsorbent (eg. silica gel) (ad)
    hads = heat of adsorption and desorption
    m1 = mass of the metallic part of the adsorber
    m2 = mass of the adsorbent
    m_cool = mass flow rate of the cooling water
    m_water = mass flow rate of the waste water
    m_ref = mass flow rate of the refrigerant
    n = adsorbent constant
    xo = maximum adsorption capacity (a constant depends on
         adsorbent-adsorbate pair)
    
    These are not part of the spec because they are determined by the control:
    tg2 = temperature at the end of the isobaric desorption
    ta2 = temperature at the end of the isobaric adsorption
    or
    x_conc = adsorption capacity at the end of the isobaric adsorption
    x_dil = adsorption capacity at the end of the isobaric desorption
    
    These are system inputs:
    t_cond = temperature of the cooling water in the condensor
    t_cool = temperature of the cooling water in the adsorber
    t_evap = temperature of the cooling water in the evaporator
    t_exhaust = temperature of the waste water
    
    These are system outputs:
    cop = coefficient of performance
    q_in = heat input by the exhaust water
    q_ref = refrigeration capacity    
    """
    def __init__(self,
        a_hex = 3.,       # m^2
        cw = 4.2,        # kJ/kg K
        c1 = 0.448,      # kJ/kg K
        c2 = 0.92,       # kJ/kg K
        cpv = 1.866,     # kJ/kg K
        hads = 2800.,     # kJ/kg
        m1 = 20.,         # kg
        m2 = 7.,            # kg
        m_water = 0.4,   #kg/s
        m_cool = 0.4,    #kg/s
        n= 0.79,
        u_des = 0.333,   # kW/m^2 k
        u_ads = 0.333,   # kW/m^2 k
        xo = 0.355,      # kg/kg
        ):
        #global a_hex  cpv cw c1 c2 end_t hads i m1 m2 m_cool m_water n start_t t_cool...
        #t_cond t_evap t_exhaust ta2 ta2_n  tg2 xo x_conc x_dil u_des u_ads
        
        self.a_hex = a_hex
        self.cw = cw
        self.c1 = c1
        self.c2 = c2
        self.cpv = cpv
        self.hads = hads
        self.m1 = m1
        self.m2 = m2
        self.m_water = m_water
        self.m_cool = m_cool
        self.n = n
        self.u_des = u_des
        self.u_ads = u_ads
        self.xo = xo

class AdsorptionChillerControl(object):
    def __init__(self,
        start_t = 0.,     # s
        end_t = 240.,     # s
        t_cool = 300.,    # K
        t_cond = 302.,    # K
        t_evap = 292.,    # K
        t_exhaust = 313., # K
        ):
        
        self.start_t = start_t
        self.end_t = end_t
        self.t_cool = t_cool
        self.t_cond = t_cond
        self.t_evap = t_evap
        self.t_exhaust = t_exhaust
        
class AdsorptionChiller(object):
    def __init__(self, spec, ctrl):
        self.spec = spec
        self.ctrl = ctrl
        self.f = Freundlich(self.spec.xo, self.spec.n)
        
    def loopOnce(self, T_d0, t):
        """Given fixed dwell time for desorption and adsorption steps,
        and temperature at start of desorption step,
        return the change in the temperature from cycle start to end.
        
        Called by solution9."""
        #global t_cond t_evap xo n
        t_cond = self.ctrl.t_cond
        t_evap = self.ctrl.t_evap
        xo = self.spec.xo
        n = self.spec.n
        
        ty, qy = self.desorption9(t, T_d0)
        T_d1 = ty[-1]
        q_d1 = qy[-1]
        P_a0 = wsr4t2p(t_evap) / ((q_d1 / xo) ** n)
        T_a0 = wsr4p2t(P_a0)
        q_a0 = q_d1
        
        tad, qad = self.adsorption9(t, T_a0)
        T_a1 = tad[-1]
        q_a1 = qad[-1]
        P_d0 = wsr4t2p(t_cond) / ((q_a1 / xo) ** n)
        T_d0_new = wsr4p2t(P_d0)
        
        dTa = T_d0_new - T_d0
        return dTa
        
    def afterSolve(self,q_low,q_high,t_cycle):
        """Called by convergeme"""
        #global m2, t_evap, t_cond, hads, c2, cw, m1, c1
        m2 = self.spec.m2
        t_evap = self.ctrl.t_evap
        t_cond = self.ctrl.t_cond
        hads = self.spec.hads
        c2 = self.spec.c2
        cw = self.spec.cw
        m1 = self.spec.m1
        c1 = self.spec.c1
        
        m_ref = 2 * (q_high - q_low) * m2 / (t_cycle)
        T_d1 = self.f.T(t_cond, q_low)
        T_a1 = self.f.T(t_evap, q_high)
        
        q_ref = m_ref * (WaterTQ2H(t_evap,1) - WaterTQ2H(t_cond,0))
        q_in = (m_ref * hads) \
            + 2 * (m2 * (c2 + q_high * cw) + m1 * c1) * (T_d1 - T_a1) / (t_cycle)
        
        if q_ref < 0:
            q_ref = 0
        
        cop = q_ref / q_in
        
        return q_ref, q_in, cop, m_ref
    
    def desorption9(self, t, T_d0):
        """Isobaric desorption process.
        
        Args:
            t : array
                A sequence of time points for which to evaluate process.
            T_d0 : float
                The initial temperature for desorption.
        
        Returns:
            ty : array
                Sequence of temperature
            qy : array
                Sequence of adsorbed mass ratio
        """
        #global tg2 n xo t_cond
        n = self.spec.n
        xo = self.spec.xo
        t_cond = self.ctrl.t_cond
        
        #[~,ty]=ode45('equation29',t, T_d0)
        # Could use scipy.integrate.quad, but that only allows one endpoint.
        # FYI scipy.integrate.ode is more generic but odeint is easy to use.
        ty = odeint(self.equation29, T_d0, t)        
        
        #qy = arrayfun(@(T)(xo*((wsr4t2p(t_cond)/wsr4t2p(T))^(1/n))), ty);
        qy = self.f.Q(t_cond, ty)
        
        #   plot(t,wsr4t2p(ty)),xlabel('Time(sec)'),ylabel('Pressure(N/m2)');
        #   hold on
        
        #tg2 = ty[-1]
        return ty, qy
        
    def desorptionFlip(self, T):
        """Integrates the equations for isobaric desorption
        over the given interval of temperature.
        
        Args:
            T : array
                Sequence of temperatures at which to evaluate process.
                Time zero corresponds to the first element of the array.
        
        Returns:
            t : array
                The time required to reach the input temperature.
            q : array
                The adsorbed mass ratio.
        """
        #global n xo t_cond
        n = self.spec.n
        xo = self.spec.xo
        t_cond = self.ctrl.t_cond
        
        #[~,t]=ode45(@(TT,tt)(1 / equation29(tt,TT)), T, 0);
        t = odeint(self.equation29flip, 0, T)
        
        #q = arrayfun(@(T)(xo*((wsr4t2p(t_cond)/wsr4t2p(T))^(1/n))), T);
        q = self.f.Q(t_cond, T)
        
        return t, q
        
    def adsorption9(self, t, T_a0):
        """Isobaric adsorption process.
        
        Args:
            t : array
                A sequence of time points for which to evaluate process.
            T_a0 : float
                The initial temperature for adsorption.
        
        Returns:
            tad : array
                Sequence of temperature
            qad : array
                Sequence of adsorbed mass ratio
        """
        #global ta2_n n xo t_evap
        n = self.spec.n
        xo = self.spec.xo
        t_evap = self.ctrl.t_evap
        
        #[~,tad]=ode45('equation49',t, T_a0);
        tad = odeint(self.equation49, T_a0, t)
        
        #qad = arrayfun(@(T)(xo*((wsr4t2p(t_evap)/wsr4t2p(T))^(1/n))), tad);
        qad = self.f.Q(t_evap, tad)
        
        # Final temperature
        ta2 = tad[-1]
        
        return tad, qad

    def adsorptionFlip(self, T):
        """Integrates the equations for isobaric adsorption
        over the given interval of temperature."""
        #global n xo t_evap
        n = self.spec.n
        xo = self.spec.xo
        t_evap = self.ctrl.t_evap
        
        #[~,t]=ode45(@(TT,tt)(1 / equation49(tt,TT)), T, 0);
        t = odeint(self.equation49flip, 0, T)
        
        #q = arrayfun(@(T)(xo*((wsr4t2p(t_evap)/wsr4t2p(T))^(1/n))), T);
        q = self.f.Q(t_evap, T)
        
        return t, q
        
    def compress(self, q):
        """Returns the dwell time (delta_t) needed to compress, assuming a
        constant adsorbate mass ratio (q).
        
        TODO: assumes U for heat transfer is u_des."""
        
        #global t_cond t_evap
        #global a_hex cw c1 c2 m1 m2 m_cool u_ads t_exhaust
        t_cond = self.ctrl.t_cond
        t_evap = self.ctrl.t_evap
        a_hex = self.spec.a_hex
        cw = self.spec.cw
        c1 = self.spec.c1
        c2 = self.spec.c2
        m1 = self.spec.m1
        m2 = self.spec.m2
        m_water = self.spec.m_water
        u_des = self.spec.u_des
        t_exhaust = self.ctrl.t_exhaust
        
        U = u_des
        T0 = self.f.T(t_evap, q)
        T1 = self.f.T(t_cond, q);
        
        NTU = (U * a_hex) / (m_water * cw);
        epsilon = 1 - exp(-NTU);
        
        A = m1 * c1 + m2 * (c2 + q * cw);
        B = m_water * cw * epsilon;
        C = t_exhaust;
        
        delta_t = A / B * log((T0 - C)/(T1 - C));
        
        return delta_t
        
    def decompress(self, q):
        """Decompression at constant adsorbed mass ratio.
        
        Returns the dwell time (delta_t) needed to decompress at the given
        adsorbed mass ratio.
        
        TODO: Assumes the U for heat transfer is u_ads."""
        
        #global t_cond t_evap
        #global a_hex cw c1 c2 m1 m2 t_cool m_cool u_ads
        t_cond = self.ctrl.t_cond
        t_evap = self.ctrl.t_evap
        a_hex = self.spec.a_hex
        cw = self.spec.cw
        c1 = self.spec.c1
        c2 = self.spec.c2
        m1 = self.spec.m1
        m2 = self.spec.m2
        m_cool = self.spec.m_cool
        u_ads = self.spec.u_ads
        t_cool = self.ctrl.t_cool
        
        U = u_ads
        T0 = self.f.T(t_cond, q);
        T1 = self.f.T(t_evap, q);
        
        NTU = (U * a_hex) / (m_cool * cw);
        epsilon = 1-exp(-NTU);
        
        A = m1 * c1 + m2 * (c2 + q * cw);
        B = m_cool * cw * epsilon;
        C = t_cool;
        
        delta_t = A / B * log((T0 - C)/(T1 - C));
        
        return delta_t
    
    def equation29(self, y, t0=None):
        """Returns the time derivative of temperature during desorption.
        
        Args:
            y : float
                Temperature in the adsorber bed. Elsewhere called ty.
            t0 : float
                Current time (ignored).
            
        Returns:
            Tdot : float
                Time derivative of temperature.
        """
        ty = y
        #global a_hex cw c1 c2 hads  m1 m2 n t_cond t_exhaust xo m_water u_des
        a_hex = self.spec.a_hex
        cw = self.spec.cw
        c1 = self.spec.c1
        c2 = self.spec.c2
        hads = self.spec.hads
        m1 = self.spec.m1
        m2 = self.spec.m2
        n = self.spec.n
        xo = self.spec.xo
        m_water = self.spec.m_water
        u_des = self.spec.u_des
        t_cond = self.ctrl.t_cond
        t_exhaust = self.ctrl.t_exhaust
        
        # TODO: extract pressure fit and mass ratio functions.
        NTU = (-u_des * a_hex) / (m_water * cw)
        epsilon = 1. - exp(-NTU)
        Qdot = epsilon * m_water * cw * (t_exhaust - ty)
        q = self.f.Q(t_cond, ty)
        term1 = m1 * c1 + m2 * (c2 + q * cw)
        term2 = m2 * hads * (xo * ((wsr4t2p(t_cond) / wsr4t2p(ty)) ** (1 / n)) \
            * ((0.7146 * ty ** 2 - 433.4 * ty + 65953) / (n * wsr4t2p(ty))))
        tdot = Qdot / (term1 + term2)
        
        return tdot
        
    def equation29flip(self, t0, y):
        """Reciprocal of equation29, to be used with odeint.
        
        Args:
            t0 : float
                The current time
            y : float
                The current adsorbed mass ratio
        
        Returns:
            The derivative of time wrt. temperature.
        """
        return 1/self.equation29(y, t0)
        
    def equation49(self, y, t0=None):
        """Returns the time derivative of temperature during adsorption.
        
        Args:
            y : float
                The temperature in the adsorption bed.
            t0 : float
                The current time (ignored)
        
        Returns:
            Tdot : float
                The temperature derivative wrt time.
        """
        tad = y
        #global a_hex cw c1 c2 cpv hads m1 m2 n t_evap t_cool xo m_cool u_ads
        a_hex = self.spec.a_hex
        cw = self.spec.cw
        c1 = self.spec.c1
        c2 = self.spec.c2
        cpv = self.spec.cpv
        hads = self.spec.hads
        m1 = self.spec.m1
        m2 = self.spec.m2
        n = self.spec.n
        xo = self.spec.xo
        m_cool = self.spec.m_cool
        u_ads = self.spec.u_ads
        t_evap = self.ctrl.t_evap
        t_cool = self.ctrl.t_cool
        
        Tdot1 = (m_cool * cw * (t_cool - tad) * (1 - exp((-u_ads * a_hex) / (m_cool * cw)))) \
            / ((m2 * (c2 + (xo * ((wsr4t2p(t_evap) / wsr4t2p(tad)) ** (1 / n)) * cw)) + m1 * c1) \
            + m2 * (hads + (cpv * (t_evap - tad))) * (xo * ((wsr4t2p(t_evap) / wsr4t2p(tad)) ** (1 / n)) * ((0.7146 * tad ** 2 - 433.4 * tad + 65953)/(n * wsr4t2p(tad)))));
        
        return Tdot1
        
    def equation49flip(self, t0, y):
        """Reciprocal of equation49 to be used with odeint.
        
        Args:
            t0 : float
                The current time.
            y : float
                The current adsorbed mass fraction.
                
        Returns:
            The derivative of time wrt. temperature."""
        return self.equation49(y, t0)
    
class Freundlich(object):
    """Store the parameters for calling Freundlich equation."""
    def __init__(self,xo,n):
        self.xo = xo
        self.n = n
    def Q(self,T_sat,T_ads):
        q = self.xo * (wsr4t2p(T_sat) / wsr4t2p(T_ads)) ** (1. / self.n)
        return q
    def T(self,T_sat,q_a):
        P_a = wsr4t2p(T_sat) / ((q_a / self.xo) ** self.n)
        T_a = wsr4p2t(P_a)
        return T_a
    def P(self, T_sat, q_a):
        P_a = wsr4t2p(T_sat) / ((q_a / self.xo) ** self.n)
        return P_a

def wsr2pt2h(p,t):
    """ WSR2PT2H  Water vapor pressure.
    h = wsr4pt2h(p,t)
        Returns enthalpy (kJ/kg) given pressure (Pa) and temperature (K).
    See also WSR4T2P."""
    h = CoolProp.CoolProp.PropsSI('H','T',t,'P',p,'Water') # J/kg
    h = h * 1e-3
    return h

def wsr4p2t(p):
    """ WSR4P2T  Water vapor saturation temperature.
    t = wsr4p2t(p) returns temperature (K) given pressure (Pa).
    
    See also WSR4P2T."""
    t = CoolProp.CoolProp.PropsSI('T','P',p,'Q',1,'Water');
    return t
    
def wsr4t2p(t):
    """WSR4T2P  Water vapor pressure.
    p = wsr4t2p(t) returns pressure (Pa) given temperature (K).
    
    See also WSR4PT2H."""
    p = CoolProp.CoolProp.PropsSI('P','T',t,'Q',1,'Water');
    return p
    
def WaterTQ2H(t,q):
    """WaterTQ2H
    h = WaterTQ2H(t,q)
        Returns enthalpy (kJ/kg) given temperature (K) and quality (kg/kg).
    See also WSR4T2P."""
    h = CoolProp.CoolProp.PropsSI('H','T',t,'Q',q,'Water'); # J/kg
    h = h * 1e-3;
    return h
    
def main():
    spec = AdsorptionChillerSpec()
    ctrl = AdsorptionChillerControl()
    chiller = AdsorptionChiller(spec, ctrl)
    t = linspace(0,1000,endpoint=True)
    T_d0 = 300
    myt = 0
    for k in range(5):
        ty,qy = chiller.desorption9(t,T_d0)
        delta_t_decompress = chiller.decompress(qy[-1])
        q_d1 = qy[-1]
        T_a0 = chiller.f.T(ctrl.t_evap, q_d1)
        tad,qad = chiller.adsorption9(t,T_a0)
        delta_t_compress = chiller.compress(qad[-1])
        plt.plot(myt + t,ty)
        myt += t[-1] + delta_t_decompress
        plt.plot(myt+t, tad)
        myt += t[-1] + delta_t_compress
        T_d0 = chiller.f.T(ctrl.t_cond, q_d1)
    plt.show()

if __name__ == "__main__":
    main()
    