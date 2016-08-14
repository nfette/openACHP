# -*- coding: utf-8 -*-
"""
Created on Mon Aug 08 07:50:51 2016

@author: nfette
"""
import CoolProp
import numpy as no
import matplotlib.pyplot as plt
from scipy.integrate import quad, odeint

class AdsorptionChillerSpec(object):
    """Variables
    a_hex = area of adsorber heat exchanger (m^2)
    cpv = specific heat of the refrigerant vapor from the evaporator
    cw = specific heat for the refrigerant( or exhaust water)
    c1 = specific heat for the metallic adsorber
    c2 = specific heat for the adsorbent(silica gel)
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
        a_hex = 3,       # m^2
        cw = 4.2,        # kJ/kg K
        c1 = 0.448,      # kJ/kg K
        c2 = 0.92,       # kJ/kg K
        cpv = 1.866,     # kJ/kg K
        hads = 2800,     # kJ/kg
        m1 = 20,         # kg
        m2=7,            # kg
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
        start_t = 0,     # s
        end_t = 240,     # s
        t_cool = 300,    # K
        t_cond = 302,    # K
        t_evap = 292,    # K
        ):
        
        self.start_t = start_t
        self.end_t = end_t
        self.t_cool = t_cool
        self.t_cond = t_cond
        self.t_evap = t_evap
        
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
        
        ty, _ = desorption9(t, T_d0)
        T_d1 = ty[-1]
        q_d1 = self.f.Q(t_cond, T_d1)
        P_a0 = wsr4t2p(t_evap) / ((q_d1 / xo) ** n)
        T_a0 = wsr4p2t(P_a0)
        q_a0 = q_d1
        
        tad, _ = adsorption9(t, T_a0)
        T_a1 = tad[-1]
        q_a1 = self.f.Q(t_evap, T_a1)
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
        
        return q_ref, q_in, cop
    
    def desorption9(self, t, T_d0):
        """Isobaric desorption"""
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
        
        tg2 = 
        return ty, qy
        
    def desorptionFlip(self, T):
        """Integrates the equations for isobaric desorption
        over the given interval of temperature."""
        global n xo t_cond
        
        [~,t]=ode45(@(TT,tt)(1 / equation29(tt,TT)), T, 0);
        q = arrayfun(@(T)(xo*((wsr4t2p(t_cond)/wsr4t2p(T))^(1/n))), T);
        return [t, q]
        
    def adsorption9(self, t, T_a0):
        """Isobaric adsorption"""
        global ta2_n n xo t_evap
        
        [~,tad]=ode45('equation49',t, T_a0);
        qad = arrayfun(@(T)(xo*((wsr4t2p(t_evap)/wsr4t2p(T))^(1/n))), tad);
        
        y3 = length(tad);
        ta2_n = real(tad(y3));
        return [tad, qad]

    def adsorptionFlip(self, T):
        """Integrates the equations for isobaric adsorption
        over the given interval of temperature."""
        global n xo t_evap
        
        [~,t]=ode45(@(TT,tt)(1 / equation49(tt,TT)), T, 0);
        q = arrayfun(@(T)(xo*((wsr4t2p(t_evap)/wsr4t2p(T))^(1/n))), T);
        return t, q
        
    def compress(self, q):
        """compress(q)
        %     Returns the dwell time needed to compress at a given adsorbate
        %     mass ratio."""
        
        global t_cond t_evap
        global a_hex cw c1 c2 m1 m2 m_cool u_ads t_exhaust
        
        T0 = Tfreundlich(t_evap, q);
        T1 = Tfreundlich(t_cond, q);
        
        NTU = (u_ads*a_hex)/(m_cool*cw);
        epsilon = 1-exp(-NTU);
        
        A = m1 * c1 + m2 * (c2 + q * cw);
        B = m_cool * cw * epsilon;
        C = t_exhaust;
        
        delta_t = A / B * log((T0 - C)/(T1 - C));
        
        return delta_t
        
    def decompress(self, q):
        """decompress(q)
            Returns the dwell time needed to decompress at a given adsorbate
            mass ratio."""
        
        global t_cond t_evap
        global a_hex cw c1 c2 m1 m2 t_cool m_cool u_ads
        
        T0 = Tfreundlich(t_cond, q);
        T1 = Tfreundlich(t_evap, q);
        
        NTU = (u_ads*a_hex)/(m_cool*cw);
        epsilon = 1-exp(-NTU);
        
        A = m1 * c1 + m2 * (c2 + q * cw);
        B = m_cool * cw * epsilon;
        C = t_cool;
        
        delta_t = A / B * log((T0 - C)/(T1 - C));
        
        return delta_t
    
    def equation29(self, _, ty):
        """Returns the time derivative of temperature during desorption."""

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
        
        tdot = (m_water * cw * (t_exhaust - ty) * (1 - exp((-u_des * a_hex) / (m_water * cw)))) / ((m2 * (c2 + (xo * ((wsr4t2p(t_cond) / wsr4t2p(ty)) ** (1 / n)) * cw)) + m1 * c1) \
            + m2 * hads * (xo * ((wsr4t2p(t_cond) / wsr4t2p(ty)) ** (1 / n)) * ((0.7146 * ty ** 2 - 433.4 * ty + 65953) / (n * wsr4t2p(ty)))));
            
        return tdot
        
    def equation49(self, _, tad):
        """Returns the time derivative of temperature during adsorption."""

        global a_hex cw c1 c2 cpv hads m1 m2 n t_evap t_cool xo m_cool u_ads
        
        Tdot1 = (m_cool*cw*(t_cool - tad)*(1-exp((-u_ads*a_hex)/(m_cool*cw))))/...
            ((m2*(c2+ (xo*((wsr4t2p(t_evap)/wsr4t2p(tad))^(1/n))*cw)) + m1*c1)+...
            m2*(hads +(cpv*(t_evap - tad)))*(xo*((wsr4t2p(t_evap)/wsr4t2p(tad))^(1/n))*((0.7146*tad^2 - 433.4*tad +65953)/(n*wsr4t2p(tad)))));
            
        return Tdot1

    
class Freundlich(object):
    """Store the parameters for calling Freundlich equation."""
    def __init__(self,xo,n):
        self.xo = xo
        self.n = n
    def Q(self,T_sat,T_ads):
        q = self.xo * (wsr4t2p(T_sat) / wsr4t2p(T_ads)) ** (1. / self.n)
        return q
    def T(self,T_sat,q_a):
        P_a = wsr4t2p(T_sat) / ((q_a / xo) ** n)
        T_a = wsr4p2t(P_a)
        return T_a

def wsr2pt2h(p,t)
    """ WSR2PT2H  Water vapor pressure.
    h = wsr4pt2h(p,t)
        Returns enthalpy (kJ/kg) given pressure (Pa) and temperature (K).
    See also WSR4T2P."""
    h = CoolProp.CoolProp.PropsSI('H','T',t,'P',p,'Water') # J/kg
    h = h * 1e-3
    return h

def wsr4p2t(p)
    """ WSR4P2T  Water vapor saturation temperature.
    t = wsr4p2t(p) returns temperature (K) given pressure (Pa).
    
    See also WSR4P2T."""
    t = CoolProp.CoolProp.PropsSI('T','P',p,'Q',1,'Water');
    return t
    
def wsr4t2p(t)
    """WSR4T2P  Water vapor pressure.
    p = wsr4t2p(t) returns pressure (Pa) given temperature (K).
    
    See also WSR4PT2H."""
    p = CoolProp.CoolProp.PropsSI('P','T',t,'Q',1,'Water');
    return p
    
def WaterTQ2H(t,q)
    """WaterTQ2H
    h = WaterTQ2H(t,q)
        Returns enthalpy (kJ/kg) given temperature (K) and quality (kg/kg).
    See also WSR4T2P."""
    h = CoolProp.CoolProp.PropsSI('H','T',t,'Q',q,'Water'); # J/kg
    h = h * 1e-3;
    return h
    
