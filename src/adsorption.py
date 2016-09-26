# -*- coding: utf-8 -*-
"""
Created on Mon Aug 08 07:50:51 2016

@author: nfette
"""
import CoolProp
#import numpy as np
from numpy import linspace, exp, log
import matplotlib.pyplot as plt
from scipy.integrate import odeint
#from hw2_1 import CelsiusToKelvin as C2K, KelvinToCelsius as K2C

class AdsorptionChillerSpec(object):
    """Variables
    
    a_hex
        [m^2] area of adsorber heat exchanger
    cpv
        [kJ/kg-K] specific heat of the refrigerant vapor from the evaporator
    cw
        [kJ/kg-K] specific heat for the refrigerant( or exhaust water)
    c1
        [kJ/kg-K] specific heat for the metallic adsorber (am)
    c2
        [kJ/kg-K] specific heat for the adsorbent (eg. silica gel) (ad)
    hads
        [kJ/kg] heat of adsorption and desorption
    m1
        [kg] mass of the metallic part of the adsorber
    m2
        [kg] mass of the adsorbent
    m_cool
        [kg/s] mass flow rate of the cooling water
    m_water
        [kg/s] mass flow rate of the waste water
    m_ref
        [kg/s] mass flow rate of the refrigerant
    n
        adsorbent constant
    xo
        [kg/kg] maximum adsorption capacity (a constant depends on
        adsorbent-adsorbate pair)
    
    These are not part of the spec because they are determined by the control:
        
    tg2
        [K] temperature at the end of the isobaric desorption
    ta2
        [K] temperature at the end of the isobaric adsorption
        
    or
    
    x_conc
        [kg/kg] adsorption capacity at the end of the isobaric adsorption
    x_dil
        [kg/kg] adsorption capacity at the end of the isobaric desorption
    """
    def __init__(self,
        a_hex = 3.,       # m^2
        cw = 4.2,        # kJ/kg K
        c1 = 0.448,      # kJ/kg K
        c2 = 0.92,       # kJ/kg K
        cpv = 1.866,     # kJ/kg K
        #cpv=0.705,
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
        
    def __repr__(self):
        import tabulate
        return tabulate.tabulate(self.__dict__.items())

class AdsorptionChillerControl(object):
    """These are system inputs:
    
    start_t 
        [s] Time at start of the cycle. Mostly useless.
    end_t (default = 240)
        [s] Time at end of the cycle.
    t_cond
        [K] temperature of the cooling water into the condensor
    t_cool
        [K] temperature of the cooling water into the adsorber
    t_evap
        [K] temperature of the cooling water into the evaporator
    t_exhaust
        [K] temperature of the waste water into the desorber"""
    
    def __init__(self,
        start_t = 0.,     # s
        end_t = 240.,     # s
        t_cool = 300.,    # K
        t_cond = 302.,    # K
        t_evap = 292.,    # K
        t_exhaust = 363., # K
        ):
        
        self.start_t = start_t
        self.end_t = end_t
        self.t_cool = t_cool
        self.t_cond = t_cond
        self.t_evap = t_evap
        self.t_exhaust = t_exhaust
    
    def __repr__(self):
        import tabulate
        return tabulate.tabulate(self.__dict__.items())
        
class AdsorptionChiller(object):
    """AdsorptionChiller provides a basic transient equilibrium model for an
    two bed adsorption cooling cycle with external heat transfer.
    
    Transient
        There are four (4) steps: desorption, decompression/expansion,
        adsorption, and compression.
    Equilibrium
        The model assumes that the adsorbent is in equilibrium with the
        refrigerant vapor in the vessel. Thus internal mass and heat transfer
        resistance are neglected, and the adsorbent is treated as lumped system
        with uniform temperature (T) and adsorbed mass ratio (q).
    Heat transfer
        Heat rates are governed by heat exchanger coefficients in the form
        of UA values, as well as stream temperature and capacity.
    
    Args
    ----
    spec
        An object inheriting AdsorptionChillerSpec that defines the internal
        parameters of the chiller
    ctrl
        An object inheriting AdsorptionChillerCtrl that defines the external
        boundary conditions and/or controls for the chiller
    
    """
    
    def __init__(self, spec, ctrl):
        self.spec = spec
        self.ctrl = ctrl
        self.f = Freundlich(self.spec.xo, self.spec.n)
        
    def loopOnce(self, T_d0, t_des, t_ads, ax1=None, ax2=None, figs=None, start_t=0):
        """Integrates through one cycle of desorption, (decompression),
        adsorption, (compression), simplified to just check convergence.
        
        Given fixed dwell time for desorption and adsorption steps,
        and temperature at start of desorption step,
        return the change in the temperature from cycle start to end.
        Optional arguments allow for plotting.

        Args
        ----
        T_d0 : float
            Temperature at start of desorption step
        t_des : array
            A sequence of time points to solve for state during integration
            of desorption
        t_ads : array
            A sequence of time points to solve for state during integration
            of adsorption
        ax1 : matplotlib.axes (optional)
            will plot the time series in the given axis
        ax2 : matplotlib.axes
            will plot the q vs T series in the given axis
        figs : matplotlib.figure
            will flush events (to watch progress of the simulation)
        start_t : float (optional)
            for the time series, specifies initial time of the cycle
        
        Returns
        -------
        dTa
            Temperature change from start of this cycle to next
            (useful to check for for convergence)
        dt
            Time to complete the entire cycle
        q_d1
            Absorbed mass ratio at end of desorption
        q_a1
            Absorbed mass ratio at end of adsorption
            
        Called by solution9."""
        #global t_cond t_evap xo n
        t_cond = self.ctrl.t_cond
        t_evap = self.ctrl.t_evap
        
        ty, qy = self.desorption9(t_des, T_d0)
        T_d1 = ty[-1]
        q_d1 = qy[-1]
        T_a0 = self.f.T(t_evap, q_d1)
        q_a0 = q_d1
        
        tad, qad = self.adsorption9(t_ads, T_a0)
        T_a1 = tad[-1]
        q_a1 = qad[-1]
        T_d0_new = self.f.T(t_cond, q_a1)
        
        dTa = T_d0_new - T_d0
        delta_t_decompress = self.decompress(q_d1)
        delta_t_compress = self.compress(q_a1)
        dt = t_des[-1] + delta_t_decompress + t_ads[-1] + delta_t_compress
        
        if ax1:
            ax1.plot(start_t + t_des, ty, 'r')
            ax1.plot(start_t + t_des[-1] + delta_t_decompress + t_ads, tad, 'b')
        if ax2:
            ax2.plot(ty,qy,'r')
            ax2.plot(tad,qad,'b')
            ax2.plot([T_d1,T_a0],[q_d1,q_a0],'o-')
            ax2.plot([T_a1,T_d0_new],[q_a1,q_a1],'o-')
        if figs:
            for fig in figs:
                try:
                    fig.canvas.draw_idle()
                    fig.canvas.flush_events()
                except:
                    fig.canvas.draw()
        
        return dTa, dt, q_d1, q_a1
        
    def afterSolve(self,q_low,q_high,t_cycle):
        """After transient simulation has converged, call this method to obtain
        time-averaged performance metrics. Assumes two beds are active.
        
        Args
        ----
        q_low
            [kg/kg] Absorbed mass ratio at end of isobaric desorption step
        q_high
            [kg/kg] Absorbed mass ratio at end of isobaric adsorption step
        t_cycle
            [s] Full cycle time for the four steps
        
        Returns
        -------
        q_ref
            [kW] cooling capacity (inbound heat flow from the cold stream)
        q_in
            [kW] inbound heat flow from the heat input stream
        cop 
            [kW/kW] thermal coefficient of performance
        m_ref
            [kg/s] mass flow rate of refrigerant
        """
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
        
        Args
        ----
        t : array
            [s] A sequence of time points for which to evaluate process.
        T_d0 : float
            [K] The initial temperature for desorption.
        
        Returns
        -------
        ty : array
            [K] Sequence of temperature
        qy : array
            [kg/kg] Sequence of adsorbed mass ratio
        """
        #global tg2 n xo t_cond
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
        
        Args
        ----
        T : array
            [K] Sequence of temperatures at which to evaluate process.
            Time zero corresponds to the first element of the array.
        
        Returns
        -------
        t : array
            [s] Sequence of the time required to reach the input temperature.
        q : array
            [kg/kg] Sequence of the adsorbed mass ratio.
        """
        #global n xo t_cond
        t_cond = self.ctrl.t_cond
        
        #[~,t]=ode45(@(TT,tt)(1 / equation29(tt,TT)), T, 0);
        t = odeint(self.equation29flip, 0, T)
        
        #q = arrayfun(@(T)(xo*((wsr4t2p(t_cond)/wsr4t2p(T))^(1/n))), T);
        q = self.f.Q(t_cond, T)
        
        return t, q
        
    def adsorption9(self, t, T_a0):
        """Isobaric adsorption process.
        
        Args
        ----
        t : array
            [s] A sequence of time points for which to evaluate process.
        T_a0 : float
            [K] The initial temperature for adsorption.
        
        Returns
        -------
        tad : array
            [K] Sequence of temperature
        qad : array
            [kg/kg] Sequence of adsorbed mass ratio
        """
        #global ta2_n n xo t_evap
        t_evap = self.ctrl.t_evap
        
        #[~,tad]=ode45('equation49',t, T_a0);
        tad = odeint(self.equation49, T_a0, t)
        
        #qad = arrayfun(@(T)(xo*((wsr4t2p(t_evap)/wsr4t2p(T))^(1/n))), tad);
        qad = self.f.Q(t_evap, tad)
        
        # Final temperature
        #ta2 = tad[-1]
        
        return tad, qad

    def adsorptionFlip(self, T):
        """Integrates the equations for isobaric adsorption
        over the given interval of temperature.
        
        Args
        ----
        T : array
            [K] Sequence of temperatures at which to evaluate process.
            Time zero corresponds to the first element of the array.
        
        Returns
        -------
        t : array
            [s] Sequence of the time required to reach the input temperature.
        q : array
            [kg/kg] Sequence of the adsorbed mass ratio.
        """
        #global n xo t_evap
        t_evap = self.ctrl.t_evap
        
        #[~,t]=ode45(@(TT,tt)(1 / equation49(tt,TT)), T, 0);
        t = odeint(self.equation49flip, 0, T)
        
        #q = arrayfun(@(T)(xo*((wsr4t2p(t_evap)/wsr4t2p(T))^(1/n))), T);
        q = self.f.Q(t_evap, T)
        
        return t, q
        
    def compress(self, q):
        """Compression of the adsorption bed at constant adsorbed mass ratio.
        
        Args
        ----
        q
            [kg/kg] Adsorbate mass ratio during the process.
        
        Returns
        -------
        delta_t
            [s] The dwell time needed to compress.
        
        TODO: assumes U for heat transfer is u_des.
        """
        
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
        T1 = self.f.T(t_cond, q)
        
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
        Assumes a heat exchanger like process.
        
        Args
        ----
        y : float
            [K] Temperature in the adsorber bed. Elsewhere called ty.
        t0 : float
            [s] Current time (ignored).
            
        Returns
        -------
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
        m_water = self.spec.m_water
        u_des = self.spec.u_des
        t_cond = self.ctrl.t_cond
        t_exhaust = self.ctrl.t_exhaust
        
        NTU = (u_des * a_hex) / (m_water * cw)
        epsilon = 1. - exp(-NTU)
        Qdot = epsilon * m_water * cw * (t_exhaust - ty)
        q = self.f.Q(t_cond, ty)
        term1 = m1 * c1 + m2 * (c2 + q * cw)
        dqdT = self.f.dQdT(t_cond, ty)
        term2 = m2 * (-hads) * dqdT
        tdot = Qdot / (term1 + term2)
        
        return tdot
        
    def equation29flip(self, t0, y):
        """Reciprocal of equation29, to be used with odeint.
        
        Args
        ----
        t0 : float
            The current time
        y : float
            The current adsorbed mass ratio
        
        Returns
        -------
        dtdT
            The derivative of time wrt. temperature.
        """
        return 1/self.equation29(y, t0)
        
    def equation49(self, y, t0=None):
        """Returns the time derivative of temperature during adsorption.
        
        Args
        ----
        y : float
            The temperature in the adsorption bed.
        t0 : float
            The current time (ignored)
        
        Returns
        -------
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
        m_cool = self.spec.m_cool
        u_ads = self.spec.u_ads
        t_evap = self.ctrl.t_evap
        t_cool = self.ctrl.t_cool
        
        NTU = (u_ads * a_hex) / (m_cool * cw)
        epsilon = (1. - exp(-NTU))
        Qdot = (m_cool * cw * (t_cool - tad) * epsilon)
        q = self.f.Q(t_evap, tad)
        dqdT = self.f.dQdT(t_evap, tad)
        term1 = m1 * c1 + m2 * (c2 + q * cw)
        term2 = m2 * ((-hads) - cpv * (t_evap - tad)) * dqdT
        Tdot1 = Qdot / (term1 + term2)
        
        return Tdot1
        
    def equation49flip(self, t0, y):
        """Reciprocal of equation49 to be used with odeint.
        
        Args
        ----
        t0 : float
            The current time.
        y : float
            The current adsorbed mass fraction.
                
        Returns
        -------
        dtdT
            The derivative of time wrt. temperature.
        """
        return 1/self.equation49(y, t0)
    
class Freundlich(object):
    """An adsorbate-refrigerant pair equilibrium model.
    Stores the parameters for calling Freundlich equation.
    
    Args
    ----
    xo
        [kg/kg] maximum adsorption capacity
    n
        exponent
    """
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
    def dQdT(self,T_sat,T_ads):
        Psat = wsr4t2p(T_ads)
        dPdT = (0.7146 * T_ads ** 2 - 433.4 * T_ads + 65953)
        q = self.Q(T_sat,T_ads)
        return -q * dPdT / (self.n * Psat)
    def dQdP(self,T_sat,T_ads):
        Psat = wsr4t2p(T_sat)
        q = self.Q(T_sat,T_ads)
        return q / (self.n * Psat)

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
    from hw2_1 import CelsiusToKelvin as C2K
    spec = AdsorptionChillerSpec()
    ctrl = AdsorptionChillerControl(t_evap=C2K(12),t_exhaust=C2K(100))
    chiller = AdsorptionChiller(spec, ctrl)
    t_ads = linspace(0,240,endpoint=True)
    t_des = linspace(0,240,endpoint=True)
    T_d0 = 300
    myt = 0
    fig1,(ax1,ax2)=plt.subplots(2,1)
    ax1.cla()
    ax2.cla()
    print(spec)
    print(ctrl)
    headers = 'Delta_T Delta_t q_low q_high'.split()
    units = 'K s kg/kg kg/kg'.split()
    fmt1 = "{:>12} "*4
    fmt2 = "{:>12.5} "*4
    print(fmt1.format(*headers))
    print(fmt1.format(*units))
    print(fmt1.format(*['----------']*4))
    for k in range(5):
        cycle = chiller.loopOnce(T_d0, t_des, t_ads, ax1, ax2, [fig1], myt)
        dT,dt,q_low,q_high = cycle
        print(fmt2.format(*cycle))
        myt += dt
        T_d0 += dT
    print(fmt1.format(*['----------']*4))
    performance = chiller.afterSolve(q_low,q_high,dt)
    
    import tabulate
    headers = 'Q_ref Q_in COP m_dot_ref'.split()
    units = 'kW kW kW/kW kg/s'.split()
    print(tabulate.tabulate(zip(headers,units,performance)))
    plt.show()
    
    return chiller,cycle,performance

if __name__ == "__main__":
    chiller,cycle,performance=main()
    