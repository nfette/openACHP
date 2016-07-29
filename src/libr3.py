# -*- coding: utf-8 -*-
"""
Created on Fri Mar 04 10:04:46 2016

@author: nfette

A single effect LiBr absorption chiller model.
"""

import CoolProp.CoolProp as CP
from hw2_1 import CelsiusToKelvin as C2K
from hw2_1 import KelvinToCelsius as K2C
import libr_props, libr_props2
import tabulate
import numpy as np
from collections import namedtuple
from scipy.optimize import fsolve
from scipy.interpolate import PchipInterpolator
import HRHX_integral_model

water = 'HEOS::Water'
librname = lambda(x): 'INCOMP::LiBr[{}]'.format(x)

ProcessPoint = namedtuple("ProcessPoint","fluid m x T P Q H D C")
def nullPP(fluid):
    return ProcessPoint(fluid,1,2,3,None,None,None,None,None)

pointType = np.dtype(dict(names="name m x T p Q h D C".split(),formats=['S32']+['d']*8))
class ProcessPoint(object):
    def __init__(self,row):
        self.row = row
def makeprop(i):
    def getx(self):
        return self.row[i]
    def setx(self,value):
        self.row[i]=value
    def delx(self):
        del self.row[i]
    return property(getx,setx,delx)
for i in pointType.names:
   setattr(ProcessPoint,i,makeprop(i))
class ProcessTable(object):
    pass
def makePointTable(names):
    table=np.zeros(len(names),dtype=pointType)
    points=ProcessTable()
    for i,name in enumerate(names):
        points.__setattr__(name,ProcessPoint(table[i]))
        table[i][0]=name
    return table,points
    
# Units in this file:
# temperature [C]
# enthalpy [J/kg]
# pressure [Pa]
# mass fraction [kg/kg total]
# effectiveness [K/K]

# The LiBr enthalpy is zero at 293.15 K, 1 atm, per
# http://www.coolprop.org/fluid_properties/Incompressibles.html#general-introduction
# We need to evaluate water enthalpy relative to that, but it is not built-in.
# http://www.coolprop.org/coolprop/HighLevelAPI.html#reference-states
T_ref = 20
P_ref = 101325
h_w_ref = CP.PropsSI('H','T',C2K(T_ref),'P',P_ref,water)

class GeneratorLiBr(object):
    """Provide process heat canonical curve for generator in various forms.
        For purpose of external heat exchange,  generator process goes ...
        P_cond, T_gen_pre, x1  (subcooled)
        P_cond, T_gen_inlet, x1 (saturated liquid) + backwards flowing vapor
        P_cond, T_gen_outlet, x2 (saturated liquid)
        
        Args:
            P : Pressure (Pa)            
                The process is assumed to occur at constant pressure level.
            m_in : Solution inlet mass flow rate (kg/s)
                The mass flow rate of the stream pumped up from low pressure.
            T_in : Solution inlet temperature (C)
                Needed to determine inlet state, which may be subcooled.
            x_in : Solution inlet mass fraction (kg/kg)
                The LiBr mass relative to total mass
            x_out : Solution outlet mass fraction (kg/kg)
                The outlet solution LiBr mass fraction
    """
        
    def __init__(self, P, m_in, T_in, x_in, x_out):        
        self.update(P, m_in, T_in, x_in, x_out)
        # Create the functions ...
        self.q = np.vectorize(self._q)
        self.T = np.vectorize(self._T)
        # Create faster functions using splines ...
        
    def update(self, P, m_in, T_in, x_in, x_out):
        self.P = P        
        self.m_in = m_in
        self.T_in = T_in
        self.x_in = x_in
        self.x_out = x_out
        
        # Find inlet enthalpy.
        # Find saturation point at inlet concentration.
        # Find ...
        self.T_sat = K2C(libr_props.temperature(self.P * 1e-5, self.x_in))
        self.T_out = K2C(libr_props.temperature(self.P * 1e-5, self.x_out))
        self.h_sat = libr_props.massSpecificEnthalpy(C2K(self.T_sat), self.x_in)
        self.h_out = libr_props.massSpecificEnthalpy(C2K(self.T_out),self.x_out)
        
        self.cp_in = libr_props.massSpecificHeat(C2K(self.T_in),self.x_in)
        # If h_in is an input:
        if False:
            preheat = (self.h_sat - self.h_in)/self.cp_in
            self.T_in = self.T_sat - preheat
        else:
            preheat = self.T_sat - self.T_in
            self.h_in = self.h_sat - preheat * self.cp_in
        
        self.h_vapor_out = CP.PropsSI('H','P', self.P,'T', C2K(self.T_sat), water) - h_w_ref
        
        # Mass balance on LiBr
        self.m_out = self.m_in * self.x_in / self.x_out
        # Mass balance on Water
        self.m_vapor_out = self.m_in - self.m_out
        self.m_total = self.m_out
        
        self.Q_preheat = self.m_in * (self.h_sat - self.h_in)
        self.Q_desorb = self.m_out * self.h_out \
                + self.m_vapor_out * self.h_vapor_out \
                - self.m_in * self.h_sat
        
        
        # Determine a reasonable limit for extending the domain.
        self.Tmax = K2C(libr_props.temperature(self.P*1e-5,libr_props.xmax))
        
        
    def __repr__(self):
        names = """P
        x_in
        x_out
        T_in
        m_in
        T_sat
        T_out
        m_out
        m_vapor_out
        h_in
        h_sat
        h_out
        h_vapor_out
        Q_preheat
        Q_desorb
        Q_total""".split()
        vals = [self.P,
                self.x_in,
                self.x_out,
                self.T_in,
                self.m_in,
                self.T_sat,
                self.T_out,
                self.m_out,
                self.m_vapor_out,
                self.h_in,
                self.h_sat,
                self.h_out,
                self.h_vapor_out,
                self.Q_preheat,
                self.Q_desorb,
                self.Q_preheat+self.Q_desorb]
        units = """Pa
        kg/kg kg/kg
        C
        kg/s
        C C
        kg/s kg/s
        J/kg J/kg J/kg J/kg
        W W W""".split()
        return tabulate.tabulate(zip(names,vals,units))

    def _q_helper(self,T):
        x_local = libr_props.massFraction(C2K(T),self.P * 1e-5)
        # Could parametrize this by x, but libr_props.temperature also has an
        # implicit solve. Only P(T,x) is explicit.
        h_vapor_local = CP.PropsSI('H',
                                   'P', self.P,
                                   'T', C2K(T), water) - h_w_ref
        h_solution_local = libr_props.massSpecificEnthalpy(C2K(T), x_local)
        # Mass balance on LiBr
        m_solution_local = self.m_in * self.x_in / x_local
        
        hlv1 = h_vapor_local - h_solution_local
        q1 = self.m_total * h_vapor_local - m_solution_local * hlv1
        return q1
        
    def _T_helper(self,q,Tguess):
        func = lambda(T):self._q_helper(T)-q
        sol = fsolve(func, Tguess)
        return sol[0]
        
    def _q(self,T):
        """Provide process heat canonical curve for generator in various forms.
        Assumes that
        * inlet solution state is obtained from self,
        * generated vapor is flowing in reverse direction and
        * it is everywhere in local equilibrium with solution.
        
        Args
        ----
        T : Temperature (deg C)
            The local temperature
        
        Returns
        -------
        q : Local progress index
            Measured as cumulative heat flux (W).
        Q : Total heat flux (W)
            The total amount of heat transferred into the generator.
        """
        
        hlv0 = self.h_vapor_out - self.h_sat
        q0 = self.m_total * self.h_vapor_out - self.m_in * hlv0
        
        if T < self.T_in:
            raise "Generator outlet must be greater than pre-inlet temperature!"
        elif T < self.T_sat:
            # The state is either saturated or subcooled.
            # Use linear interpolation.
            q = self.m_in * self.cp_in * (T - self.T_in)
        else:
            x_local = libr_props.massFraction(C2K(T),self.P * 1e-5)
            h_vapor_local = CP.PropsSI('H',
                                       'P', self.P,
                                       'T', C2K(T), water) - h_w_ref
            h_solution_local = libr_props.massSpecificEnthalpy(C2K(T), x_local)
            # Mass balance on LiBr
            m_solution_local = self.m_in * self.x_in / x_local
            
            hlv1 = h_vapor_local - h_solution_local
            q1 = self.m_total * h_vapor_local - m_solution_local * hlv1
            
            q = (q1 - q0) + self.Q_preheat
        # TODO
        # If q > Q (iff T > T_solution_outlet) then result is invalid.
        return q
        
    def _T(self,q):
        if q < 0:
            raise "Process index must be greater than 0!"
        elif q < self.Q_preheat:
            T = self.T_in + q / (self.m_in * self.cp_in)
        #elif q < self.Q_desorb + self.Q_preheat:
        else:
            hlv0 = self.h_vapor_out - self.h_sat
            q0 = self.m_total * self.h_vapor_out - self.m_in * hlv0
            q1 = q0 + (q - self.Q_preheat)
            alpha = (q - self.Q_preheat) / self.Q_desorb
            Tguess = alpha*self.T_sat+(1-alpha)*self.T_out
            T = self._T_helper(q1,Tguess)
        #else:
        #    raise "Process index extends past process
        return T
        
class GeneratorLiBrInterpolated(GeneratorLiBr):
    def __init__(self, P, m_in, T_in, x_in, x_out, debug=False):
        print "I am still alive!"
        self.update(P, m_in, T_in, x_in, x_out)
        TT = np.linspace(self.T_in,self.Tmax)
        qfunc = np.vectorize(self._q)
        qq = qfunc(TT)
        if np.isnan(qq).any():
            print "There is a problem with some nans"

        # Create extended domain such that interpolate saturates at endpoints.        
        qq1 = np.resize(qq,qq.size+1)
        TT1 = np.resize(TT,TT.size+1)
        qq1[-1] = qq1[-2]
        TT1[-1] = TT1[-2] + 1
        if (np.diff(qq1) < 0).any():
            print "Captain, it's a non-monotonic function!"        
        self.q = PchipInterpolator(TT1,qq1,extrapolate=True)
        
        # Need to use fresh arrays because they are referenced.
        qq2 = np.resize(qq,qq.size+1)
        TT2 = np.resize(TT,TT.size+1)
        qq2[-1] = qq2[-2] * 1.02
        TT2[-1] = TT2[-2]
        self.T = PchipInterpolator(qq2,TT2,extrapolate=True)

        # Show that it worked
        if debug:
            print tabulate.tabulate(zip(TT1,qq1))
            print tabulate.tabulate(zip(TT2,qq2))
            import matplotlib.pyplot as plt
            plt.figure()                
            plt.plot(TT1,qq1,'.'); plt.title("qqmod vs TTmod for q()")
            plt.figure()
            plt.plot(qq2,TT2,'.'); plt.title("TTmod vs qqmod for T()")
        
class AbsorberLiBr1(object):
    """Provides a canonical heat (output) curve for a LiBr water vapor absorber.
    
    Assumes:
        Vapor comes into equilibrium with surface only where it is absorbed.
    
    Inputs
    ------
        P : number (Pa)
            Process pressure
        m_in : number (kg/s)
            Solution inlet mass flow rate
        h_in : number (J/kg)
            Solution inlet enthalpy
        x_in : number (kg/kg)
            Solution inlet mass fraction LiBr
        T_vapor_inlet : number (deg C)
            Vapor inlet temperature
    """
    def __init__(self,P,m_in,h_in,x_in,h_vapor_inlet,debug=False):
        self.P=P
        self.m_in=m_in
        self.h_in=h_in
        self.x_in=x_in
        self.h_vapor_inlet = h_vapor_inlet
        
        self.T_sat = K2C(libr_props.temperature(self.P * 1e-5,
                                                          self.x_in))
        self.h_sat = libr_props.massSpecificEnthalpy(C2K(self.T_sat),self.x_in)
        
        if self.h_in > self.h_sat:
            q,t,xl = libr_props.twoPhaseProps(self.h_in,self.P*1e-5,self.x_in)
            # Pre-cooling is required to reach saturation temperature
            self.Q_pre_cool = self.m_in * (self.h_sat - self.h_in)
            self.T_in = K2C(t)
            prepoints_h = np.linspace(self.h_in, self.h_sat)
            prepoints_T = np.zeros_like(prepoints_h)
            prepoints_q = np.zeros_like(prepoints_h)
            prepoints_x = np.zeros_like(prepoints_h)
            for i,h in enumerate(prepoints_h):
                q,t,xl = libr_props.twoPhaseProps(h,self.P*1e-5,self.x_in)
                prepoints_T[i] = K2C(t)
                prepoints_q[i] = self.m_in * (h - self.h_in)
                prepoints_x[i] = xl
        else:
            self.Q_pre_cool = 0
            self.T_in = K2C(libr_props.temperature(self.P * 1e-5, self.x_in))
            prepoints_T = []
            prepoints_q = []
            prepoints_x = []
        
        # Set up bounds and points for interpolation.
        # Absorber limit is liquid water.
        self.Tmin = CP.PropsSI('T','P',P,'Q',0,water)
        x_points = np.linspace(x_in,0.1)
        T_points = np.zeros_like(x_points)
        q_points = np.zeros_like(x_points)
        #q_func = np.vectorize(self._q)
        #q_points = q_func(T_points)
        for i,x in enumerate(x_points):
            T_points[i],q_points[i] = self._qx(x)

        x_points = np.concatenate([prepoints_x,x_points])
        T_points = np.concatenate([prepoints_T,T_points])
        q_points = np.concatenate([prepoints_q,q_points])
        
        # Forward function
        T_points1 = np.resize(T_points,len(T_points)+1)
        q_points1 = np.resize(q_points,len(T_points)+1)
        T_points1[-1] = T_points[-1] - 5
        q_points1[-1] = q_points[-1]
        
        # Inverse function
        T_points2 = np.resize(T_points,len(T_points)+1)
        q_points2 = np.resize(q_points,len(T_points)+1)
        T_points2[-1] = T_points[-1]
        q_points2[-1] = q_points[-1] * 1.05
        
        if debug:        
            import matplotlib.pyplot as plt
            import tabulate
            xmod = np.resize(x_points,len(x_points)+1)
            print tabulate.tabulate(zip(xmod,T_points1,q_points1,
                                    T_points2,q_points2),
                                    headers=['x','T1','q1','T2','q2'])
            plt.figure(); plt.plot(T_points1,q_points1); plt.title("q(T)")
            plt.figure(); plt.plot(q_points2,T_points2); plt.title("T(q)")
        # Interpolate data must be in increasing order, so we reverse it
        # compared to reaction direction.
        self.q = PchipInterpolator(T_points1[::-1], q_points1[::-1])
        self.T = PchipInterpolator(q_points2[::-1], T_points2[::-1])
        
    def _q(self,T):
        # First, determine the mass fraction here.
        if T > self.T_sat:
            raise ValueError("This function is for subcooled LiBr...")
        else:
            x_local = libr_props.massFraction(C2K(T),self.P*1e-5)
        return self._qx(x_local)
    
    def _qx(self,x_local):
        T = K2C(libr_props.temperature(self.P*1e-5,x_local))
        dx = self.x_in - x_local
        # TODO
        h_local = libr_props.massSpecificEnthalpy(C2K(T),x_local)
        # And calculate
        m_vapor = self.m_in * dx / x_local
        m_out = self.m_in + m_vapor
        # Heat is inbound.
        result = -self.m_in * self.h_in - m_vapor * self.h_vapor_inlet \
            + m_out * h_local + self.Q_pre_cool
        return T,result
        
    def __repr__(self):
        names = """P
        m_in
        x_in
        T_in
        T_vapor_inlet
        h_in
        h_vapor_inlet""".split()
        vals = [self.P,
                self.m_in,
                self.x_in,
                self.T_in,
                self.T_vapor_inlet,
                self.h_in,
                self.h_vapor_inlet]
        units = """Pa
        kg/kg kg/kg
        C
        kg/s
        C C
        kg/s kg/s
        J/kg J/kg J/kg J/kg
        W W W""".split()
        return tabulate.tabulate(zip(names,vals,units))
        
class ChillerLiBr1(object):
    def __init__(self,
                 T_evap=1.5, T_cond=39.9,
                 x1=0.567, x2=0.624,
                 Eff_SHX=0.64, m_pump=0.05):
        """Args
        ----
            T_evap : float
                Evaporator saturation temperature (deg C)
            T_cond : float
                Condenser saturation temperature (deg C)
            x1 : float
                Pump side mass fraction of LiBr in stream (low) [kg/kg]
            x2 : float
                Return side mass fraction of LiBr in stream (high/concentrate)
            Eff_SHX : float
                Effectiveness of the solution heat exchanger (K/K), 0 to 1
            m_pump : float
                Mass flow rate through the solution pump (kg/s)
        """
        self.T_evap = T_evap
        self.T_cond = T_cond
        self.x1 = x1
        self.x2 = x2
        self.m_pump = m_pump
        self.Eff_SHX = Eff_SHX
        self.dx = x1 - x2
        
        self.P_evap = CP.PropsSI('P','T',C2K(T_evap),'Q',1,water)
        self.P_cond = CP.PropsSI('P','T',C2K(T_cond),'Q',1,water)
        
        self.stateLabels = """abs_outlet
pump_outlet
gen_inlet
gen_sat_liquid
gen_outlet
SHX_conc_outlet
abs_inlet
abs_sat_liquid
gen_vapor_outlet
cond_sat_vapor
cond_outlet
evap_inlet
evap_sat_liquid
evap_sat_vapor
evap_outlet""".split('\n')
        #self.states=dict((k, nullPP('LiBrH2O')) for k in self.stateLabels)
        self.stateTable,self.states=makePointTable(self.stateLabels)
        
        
        self.T_gen_inlet = 0
        self.T_gen_outlet = 0
        self.T_abs_inlet_max = 0
        self.T_abs_outlet_max = 0
        self.h_gen_inlet = 0
        self.h_gen_outlet = 0
        self.h_abs_inlet = 0
        self.h_abs_outlet = 0
        self.m_concentrate = 0
        self.m_refrig = 0
        self.T_SHX_concentrate_outlet = 0
        self.Q_SHX = 0
        self.T_abs_pre = np.nan
        self.h_abs_pre = np.nan
        self.Q_abs_pre_cool = 0
        self.P_abs_pre = np.nan
        self.Q_abs_main = 0
        self.Q_abs_total = 0
        self.T_gen_pre = np.nan
        self.Q_gen_pre_heat = 0
        self.Q_gen_main = 0
        self.Q_gen_total = 0
        self.Q_condenser_reject = 0
        self.Q_evap_heat = 0
        self.COP = 0
        self.W_pump = 0
        self.f = np.inf
        self.x_abs_pre = self.x2
    
    # These routines allow updating solution
    def setT_evap(self,T_evap):
        self.T_evap = T_evap
        self.P_evap = CP.PropsSI('P','T',C2K(T_evap),'Q',1,water)
    def setT_cond(self,T_cond):
        self.T_cond = T_cond
        self.P_cond = CP.PropsSI('P','T',C2K(T_cond),'Q',1,water)
        
    def ZeroCheck(self):
        return self.W_pump + self.Q_evap_heat + self.Q_gen_total - self.Q_condenser_reject - self.Q_abs_total
        
    def iterate1(self):
        """Update the internal parameters."""
        self.T_gen_inlet = K2C(libr_props.temperature(self.P_cond*1e-5,
                                                      self.x1))
        self.T_gen_outlet = K2C(libr_props.temperature(self.P_cond * 1e-5,
                                                   self.x2))
        self.T_abs_inlet_max = K2C(libr_props.temperature(self.P_evap * 1e-5,
                                                          self.x2))
        self.T_abs_outlet_max = K2C(libr_props.temperature(self.P_evap * 1e-5,
                                                           self.x1))
        
        self.h_gen_inlet = libr_props.massSpecificEnthalpy(
            C2K(self.T_gen_inlet), self.x1)
        self.h_gen_outlet = libr_props.massSpecificEnthalpy(
            C2K(self.T_gen_outlet), self.x2)
        self.h_abs_inlet = libr_props.massSpecificEnthalpy(
            C2K(self.T_abs_inlet_max), self.x2)
        self.h_abs_outlet = libr_props.massSpecificEnthalpy(
            C2K(self.T_abs_outlet_max), self.x1)
        
        # Mass balance on LiBr
        self.m_concentrate = self.m_pump * self.x1 / self.x2
        # Mass balance on Water
        self.m_refrig = self.m_pump - self.m_concentrate
        self.f = self.m_pump / self.m_refrig

        # Compute SHX outlets, assuming concentrate limits heat flow (C_min)
        # Neglect pump work for the present.
        DeltaT_max = self.T_gen_outlet - self.T_abs_outlet_max
        DeltaT_SHX_concentrate = self.Eff_SHX * DeltaT_max
        self.T_SHX_concentrate_outlet = self.T_gen_outlet \
            - DeltaT_SHX_concentrate
        self.h_SHX_concentrate_outlet = libr_props.massSpecificEnthalpy(
            C2K(self.T_SHX_concentrate_outlet), self.x2)
        self.Q_SHX = self.m_concentrate \
            * (self.h_gen_outlet - self.h_SHX_concentrate_outlet)
        
        # Expansion valve
        self.h_abs_pre = self.h_SHX_concentrate_outlet
        if self.h_abs_pre > self.h_abs_inlet:
            # Pre-cooling is required to reach saturation temperature
            self.Q_abs_pre_cool = self.m_concentrate \
                * (self.h_abs_pre - self.h_abs_inlet)
            q,t,xl = libr_props.twoPhaseProps(self.h_abs_pre,
                                              self.P_evap*1e-5,
                                              self.x2)
            self.T_abs_pre = K2C(t)
            self.x_abs_pre = xl
            # ignore vapor quality, q
            # Minimum vapor pressure for absorption to occur
            self.P_abs_pre = np.inf
        else:
            self.Q_abs_pre_cool = 0
            #self.T_abs_pre = K2C(CP.PropsSI('T',
            #    'H', self.h_abs_pre,
            #    'P', self.P_evap,
            #    librname(self.x2)))
            self.T_abs_pre = np.nan
            # Minimum vapor pressure for absorption to occur
#            self.P_abs_pre = CP.PropsSI('P',
#                'T', C2K(self.T_abs_pre),
#                'Q', 0,
#                librname(self.x2))
            self.P_abs_pre = np.nan
                
        # Heat rejection in absorber: energy balance
        self.h_abs_vapor_inlet = CP.PropsSI('H',
            'P',self.P_evap,
            'Q',1,
            water) - h_w_ref
        self.Q_abs_main = self.m_refrig * self.h_abs_vapor_inlet \
            + self.m_concentrate * self.h_abs_inlet \
            - self.m_pump * self.h_abs_outlet
        self.Q_abs_total = self.Q_abs_main + self.Q_abs_pre_cool
        
        # Energy balance in SHX, pump side
        D_in = CP.PropsSI('D',
            'T',C2K(self.T_abs_outlet_max),
            'Q',0,
            librname(self.x1))
        DeltaH_pump = (self.P_cond - self.P_evap) / D_in
        self.W_pump = self.m_pump * DeltaH_pump
        self.h_pump_outlet = self.h_abs_outlet + DeltaH_pump
        DeltaH_SHX_pumpside = self.Q_SHX / self.m_pump
        self.h_gen_pre = self.h_pump_outlet + DeltaH_SHX_pumpside
        if self.h_gen_pre > self.h_gen_inlet:
            # Flash steam
            self.T_gen_pre = np.nan
        else:
            # The state is either saturated or subcooled.
            # We need to calculate the temperature from specific heat.
            cp = libr_props.massSpecificHeat(C2K(self.T_gen_inlet),self.x1)
            deltaHsub = self.h_gen_inlet - self.h_gen_pre
            deltaT = deltaHsub / cp
            self.T_gen_pre = self.T_gen_inlet - deltaT
        
        self.Q_gen_pre_heat = self.m_pump * (self.h_gen_inlet - self.h_gen_pre)
        
        # Heat input to generator: energy balance
        self.h_gen_vapor_outlet = CP.PropsSI('H',
            'P', self.P_cond,
            'T', C2K(self.T_gen_inlet), water) - h_w_ref
        self.vapor_superheat = self.T_gen_inlet - self.T_cond
        self.Q_gen_main = self.m_refrig * self.h_gen_vapor_outlet \
            + self.m_concentrate * self.h_gen_outlet \
            - self.m_pump * self.h_gen_inlet
        self.Q_gen_total = self.Q_gen_main + self.Q_gen_pre_heat
        
        # Condenser
        self.h_condenser_outlet = CP.PropsSI('H',
            'P',self.P_cond,
            'Q', 0, water) - h_w_ref
        self.Q_condenser_reject = self.m_refrig * (self.h_gen_vapor_outlet
            - self.h_condenser_outlet)
        
        # Expansion valve
        self.h_evap_inlet = self.h_condenser_outlet
        
        # Evaporator
        self.h_evap_outlet = self.h_abs_vapor_inlet
        self.Q_evap_heat = self.m_refrig * (self.h_evap_outlet
            - self.h_evap_inlet)
        
        self.COP = self.Q_evap_heat / self.Q_gen_total
    
    def updateGenerator(self,Q_gen):
        genStream = self.getGeneratorStream()
        self.h_gen_inlet = genStream.h_sat
        self.T_gen_inlet = genStream.T_sat
        
        if Q_gen < genStream.Q_preheat:
            self.T_gen_outlet = genStream.T(Q_gen)
            self.x2 = self.x1
            # error?
        else:
            self.T_gen_outlet = genStream.T(Q_gen)
            self.x2 = libr_props.massFraction(C2K(self.T_gen_outlet),
                                              self.P_cond)
            genStream.update(self.P_cond, self.m_pump, self.T_gen_pre,
                             self.x1, self.x2)
            self.h_gen_outlet = genStream.h_out
            self.h_gen_vapor_outlet = genStream.h_vapor_out
        self.dx = self.x1 - self.x2
        # Mass balance on LiBr
        self.m_concentrate = self.m_pump * self.x1 / self.x2
        # Mass balance on Water
        self.m_refrig = self.m_pump - self.m_concentrate
        self.f = self.m_pump / self.m_refrig
            
    def updateSHX_hot_side(self):
        DeltaT_max = self.T_gen_outlet - self.T_abs_outlet_max
        DeltaT_SHX_concentrate = self.Eff_SHX * DeltaT_max
        self.T_SHX_concentrate_outlet = self.T_gen_outlet \
            - DeltaT_SHX_concentrate
        self.h_SHX_concentrate_outlet = libr_props.massSpecificEnthalpy(
            C2K(self.T_SHX_concentrate_outlet), self.x2)
        self.Q_SHX = self.m_concentrate \
            * (self.h_gen_outlet - self.h_SHX_concentrate_outlet)
        
        
    def getHeatCurve(self):
        """Returns (Heat,T), arrays showing progress of the process.
        Note: absorber heat input is negative.
        Learn this from (my revision to) example 3.3.
        """
        Q, result = 0, []
        # Starting coordinates
        result.append((0,self.T_abs_pre))
        # Heat input to reach a saturated state.
        Q += self.m_concentrate * (self.h_abs_inlet - self.h_abs_pre)
        result.append((Q,self.T_abs_inlet_max))
        # Heat input to reach outlet.
        Q += self.m_pump * self.h_abs_outlet \
            - self.m_concentrate * self.h_abs_inlet \
            - self.m_refrig * self.h_abs_vapor_inlet
        result.append((Q,self.T_abs_outlet_max))
        # Pump -- no heat, just work
        
        # Pump outlet to generator through SHX
        Q += self.m_pump * self.h_gen_pre - self.m_pump * self.h_abs_outlet
        result.append((Q,self.T_gen_pre))
        # Generator preheat
        Q += self.m_concentrate * self.h_gen_inlet \
                - self.m_concentrate * self.h_gen_pre
        result.append((Q,self.T_gen_inlet))        
        # Generator proper
        Q += self.m_concentrate * self.h_gen_outlet \
                + self.m_refrig * self.h_gen_vapor_outlet \
                - self.m_pump * self.h_gen_inlet
        result.append((Q,self.T_gen_outlet))
        # SHX, concentrate side
        Q += self.m_concentrate * self.h_SHX_concentrate_outlet \
                - self.m_concentrate * self.h_gen_outlet
        result.append((Q,self.T_SHX_concentrate_outlet))
        # Solution expander
        result.append((Q,self.T_abs_pre))
        # Condenser cool to saturated
        result.append((Q,self.T_gen_inlet))
        h_condenser_sat = CP.PropsSI("H","T",C2K(self.T_cond),"Q",0,
                                     "HEOS::water")
        Q += self.m_refrig * h_condenser_sat \
                - self.m_refrig * self.h_gen_vapor_outlet
        result.append((Q,self.T_cond))
        # Real condense
        Q += self.m_refrig * (self.h_condenser_outlet - h_condenser_sat)
        result.append((Q,self.T_cond))
        # What if condenser subcools? Later.
        # Expander
        T_into_evap = K2C(CP.PropsSI("T","H",self.h_condenser_outlet,
                                 "P",self.P_evap,"HEOS::water"))
        result.append((Q,T_into_evap))
        # Evaporator
        Q += self.m_refrig * (self.h_evap_outlet - self.h_evap_inlet)
        result.append((Q,self.T_evap))
        
        return zip(*result)

    def iterate2(self,T_gen_outlet,T_abs_outlet):
        """Resolve the concentrations. Not yet implemented."""
        pass
    
    def buildGeneratorHeatCurve(self):
        """Provide process heat canonical curve for generator in various forms.
        For purpose of external heat exchange,  generator process goes ...
        P_cond, T_gen_pre, x1  (subcooled)
        P_cond, T_gen_inlet, x1 (saturated liquid) + backwards flowing vapor
        P_cond, T_gen_outlet, x2 (saturated liquid)
        
        Returns
        -------
        T : callable
            Maps process progress in terms of heat flux (W) to local
            temperature T (deg C).
        q : callable
            Maps process local temperature T (deg C) and total heat flux Q (W)
            to the local progress index, cumulative heat flux q (W).
        """
        self.genpoints, Q = [], 0
        # 0, pre-inlet
        self.genpoints.append((Q,self.T_gen_pre))
        # 1, Generator preheat
        Q += self.m_concentrate * self.h_gen_inlet \
                - self.m_concentrate * self.h_gen_pre
        self.genpoints.append((Q,self.T_gen_inlet))        
        # 2, Generator proper
        Q += self.m_concentrate * self.h_gen_outlet \
                + self.m_refrig * self.h_gen_vapor_outlet \
                - self.m_pump * self.h_gen_inlet
        self.genpoints.append((Q,self.T_gen_outlet))
        
    def generatorHeatCurveQ(self,T,x_out):
        """Provide process heat canonical curve for generator in various forms.
        Assumes that
        * inlet solution state is obtained from self,
        * generated vapor is flowing in reverse direction and
        * it is everywhere in local equilibrium with solution.
        
        Args
        ----
        T : Temperature (deg C)
            The local temperature
        x_out : Mass fraction (kg/kg)
            The outlet solution LiBr mass fraction
        
        Returns
        -------
        q : Local progress index
            Measured as cumulative heat flux (W).
        Q : Total heat flux (W)
            The total amount of heat transferred into the generator.
        """
        
        # We may need m_concentrate. It depends on Q -- solve from inlet.
        h_solution_inlet = self.h_gen_inlet
        h_vapor_outlet = self.h_gen_vapor_outlet
        T_solution_outlet = K2C(libr_props.temperature(self.P_cond * 1e-5,
                                               x_out))
        h_solution_outlet = libr_props.massSpecificEnthalpy(
            C2K(T_solution_outlet), x_out)
        # Mass balance on LiBr
        m_solution_outlet = self.m_pump * self.x1 / x_out
        # Mass balance on Water
        m_vapor_outlet = self.m_pump - m_solution_outlet
        m_total = m_solution_outlet
        
        hlv0 = h_vapor_outlet - h_solution_inlet
        q0 = m_total * h_vapor_outlet - self.m_pump * hlv0
        
        Q = m_solution_outlet * h_solution_outlet \
            + m_vapor_outlet * h_vapor_outlet \
            - self.m_pump * h_solution_inlet
        
        q = 0
        T0,T1,T2=self.genpoints[0][1],self.genpoints[1][1],self.genpoints[2][1]
        Q0,Q1,Q2=self.genpoints[0][0],self.genpoints[1][0],self.genpoints[2][0]
        if T < T0:
            raise "Generator outlet must be greater than pre-inlet temperature!"
        elif T < T1:
            # The state is either saturated or subcooled.
            # Use linear interpolation.
            q = (T - T0) / (T1 - T0) * (Q1 - Q0)
        else:
            x_local = libr_props.massFraction(C2K(T),self.P_cond * 1e-5)
            h_vapor_local = CP.PropsSI('H',
                                       'P', self.P_cond,
                                       'T', C2K(T), water) - h_w_ref
            h_solution_local = libr_props.massSpecificEnthalpy(C2K(T), x_local)
            # Mass balance on LiBr
            m_solution_local = self.m_pump * self.x1 / x_local
            
            hlv1 = h_vapor_local - h_solution_local
            q1 = m_total * h_vapor_local - m_solution_local * hlv1
            
            q = (q1 - q0) + (Q1 - Q0)
        # TODO
        # If q > Q (iff T > T_solution_outlet) then result is invalid.
        return q, Q
        
    def generatorHeatCurveT(self,Q):
        """Provide process heat canonical curve for generator in various forms."""
        T = 0
        return T
    
    def getGeneratorStream(self):
        gen = GeneratorLiBrInterpolated(self.P_cond,self.m_pump,self.T_gen_pre,
                                        self.x1,self.x2)
        return gen
        
    def getAbsorberStream(self):
        absorber = AbsorberLiBr1(self.P_evap,
                                 self.m_concentrate,
                                 self.h_abs_pre,
                                 self.x_abs_pre,
                                 self.h_evap_outlet)
        return absorber
    
    def getCondenserStream(self):
        h_rel = self.h_gen_vapor_outlet + h_w_ref
        condenser = HRHX_integral_model.waterStream(self.P_cond,h_rel,self.m_refrig)
        return condenser
        
    def getEvaporatorStream(self):
        h_rel = self.h_evap_inlet + h_w_ref
        evaporator = HRHX_integral_model.waterStream(self.P_evap,h_rel,self.m_refrig)
        return evaporator
    
    def __repr__(self):
        names = """T_evap T_cond P_evap P_cond
        x1 x2
        T_gen_inlet T_gen_outlet T_abs_inlet_max T_abs_outlet_max
        h_gen_inlet h_gen_outlet h_abs_inlet h_abs_outlet
        m_pump m_concentrate m_refrig
        Eff_SHX
        T_SHX_concentrate_outlet Q_SHX
        T_abs_pre h_abs_pre x_abs_pre Q_abs_pre_cool P_abs_pre
        Q_abs_main Q_abs_total
        T_gen_pre
        Q_gen_pre_heat Q_gen_main Q_gen_total
        Q_condenser_reject Q_evap_heat COP
        W_pump
        self.f
        ZeroCheck
        string""".split()
        vals = [self.T_evap, self.T_cond, self.P_evap, self.P_cond,
            self.x1, self.x2,
            self.T_gen_inlet, self.T_gen_outlet, self.T_abs_inlet_max, self.T_abs_outlet_max,
            self.h_gen_inlet, self.h_gen_outlet, self.h_abs_inlet, self.h_abs_outlet,
            self.m_pump, self.m_concentrate, self.m_refrig,
            self.Eff_SHX,
            self.T_SHX_concentrate_outlet, self.Q_SHX,
            self.T_abs_pre, self.h_abs_pre, self.x_abs_pre,
            self.Q_abs_pre_cool, self.P_abs_pre,
            self.Q_abs_main, self.Q_abs_total,
            self.T_gen_pre,
            self.Q_gen_pre_heat, self.Q_gen_main, self.Q_gen_total,
            self.Q_condenser_reject, self.Q_evap_heat, self.COP,
            self.W_pump,
            self.f,
            self.ZeroCheck(),
            False]
        units = """C C Pa Pa
        kg/kg kg/kg
        C C C C
        J/kg J/kg J/kg J/kg
        kg/s kg/s kg/s
        K/K
        C W
        C J/kg kg/kg W Pa
        W W
        C
        W W W
        W W W/W
        W
        kg/kg
        W
        none""".split()
        vartable = tabulate.tabulate(zip(names,vals,units))
        statetable = tabulate.tabulate(self.stateTable,pointType.names)
        return vartable + "\n" + statetable
        
def main():
    if True:
        # Example 6.1 in the book
        P1,P2 = 673, 7445
        T1 = K2C(CP.PropsSI('T','P',P1,'Q',1,water))
        T2 = K2C(CP.PropsSI('T','P',P2,'Q',1,water))
        # Trying different inputs
        T1, T2 = 1, 35
        c = ChillerLiBr1(T1,T2,0.5,0.7)
        c.x2=libr_props2.Xsat(89.9,c.P_cond)
        c.x1=libr_props2.Xsat(32.7,c.P_evap)
        # Custom example
        c = ChillerLiBr1(T_evap=5,T_cond=45,x1=0.6026,x2=0.66)
        print "Initializing..."
        print c
        print "Iterating..."
        try:
            c.iterate1()
        finally:
            print c
        
    if True:
        # Figure 6.3 in the book
        Eff_SHX = np.linspace(0,1)
        COP = np.zeros_like(Eff_SHX)
        for i in range(len(Eff_SHX)):
            c = ChillerLiBr1(Eff_SHX=Eff_SHX[i])
            try:
                c.iterate1()
                COP[i] = c.COP
            except:
                pass
        if False:
            import matplotlib.pyplot as plt
            plt.plot(Eff_SHX, COP)
            plt.show()
    
    return c
    
if __name__ == "__main__":
    c=main()