# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 13:23:47 2016

A single effect ammonia-water cycle.

@author: nfette
"""

import CoolProp.CoolProp as CP
from hw2_1 import CelsiusToKelvin as C2K
from hw2_1 import KelvinToCelsius as K2C
import tabulate
import numpy as np
import scipy.interpolate
import scipy.optimize
import HRHX_integral_model

from ammonia_props import AmmoniaProps, StateType, convert_state_list_to_array, CStateTable

amm = AmmoniaProps()

class stateIterator:
    def __init__(self, chiller):
        self.chiller = chiller
        self.i = iter(chiller.points)

    def __iter__(self):
        return self

    def __next__(self):
        return self.chiller.__getattribute__(self.i.__next__())


class CVariablesTable:

    def __init__(self, names, values, units=None):
        """A convenience class for a list of variable names, units, and values.
        TODO: utilize existing libraries ... maybe pandas.
        """
        dt = [('name', 'S256'), ('unit', 'S256'), ('value', 'f')]
        self.table = np.zeros_like(names, dtype=dt)
        self.table['name'] = names
        self.table['value'] = values
        self.table['unit'] = units

    def tabulate(self, **kwargs):
        return tabulate.tabulate(self.table,
             self.table.dtype.names,
             **kwargs
        )

    def __repr__(self):
        return self.tabulate()

    def _repr_html_(self):
        return self.tabulate(tablefmt="html")


class AmmoniaAbsorberStream(object):
    """TODO: Account for non-saturated solution inlet."""

    def __init__(self, weak_inlet, refrig_inlet, m_weak, debug=False):
        self.weak_inlet = weak_inlet
        self.refrig_inlet = refrig_inlet
        self.m_weak = m_weak

        T_points = []
        q_points = []
        x_points = []
        self.q_pre = 0

        if weak_inlet.isSubcooled():
            #print("Note: Absorber inlet is subcooled")

            # If so, then any absorption prior to saturation does not require
            # heat transfer to proceed. By design, this should probably not
            # happen; instead, expand to lower pressure until saturation is
            # reached. Anyway, we should calculate the saturated state where
            # heat transfer starts. Assume:
            #   1. Pre-saturation absorption is at constant temperature.
            self.sat_inlet = amm.props2(T=weak_inlet.T, P=weak_inlet.P, Qu=0)
            x_sat = self.sat_inlet.x
            # m_sat - m_refrig = m_weak
            # x_sat m_sat - m_refrig x_refrig = m_weak x_weak
            self.m_sat = m_weak * (weak_inlet.x - refrig_inlet.x) / (x_sat - refrig_inlet.x)
            m_refrig = m_weak * (weak_inlet.x - x_sat) / (x_sat - refrig_inlet.x)

        elif weak_inlet.Qu > 0:
            #print("Note: Absorber inlet contains some vapor.")

            # Inlet is already a (super)saturated state. However, there may be
            # some vapor at the inlet. So, we should cool the thing down.
            self.sat_inlet = amm.props2(P=weak_inlet.P, Qu=0, x=weak_inlet.x)
            self.m_sat = m_weak
            if weak_inlet.isSuperheated():
                #print("Note: Absorber inlet is superheated.")

                # If so, then the fluid must be further cooled to saturation
                # temperature before any vapor can be absorbed. This should not
                # have happened, by design! Assume:
                #   1. We must cool to saturated vapor.
                vapor = amm.props2(P=weak_inlet.P, Qu=1, x=weak_inlet.x)

                prepoints_h = np.linspace(weak_inlet.h, vapor.h, 10, endpoint=False)
                for i, h in enumerate(prepoints_h):
                    state = amm.props2(P=weak_inlet.P, x=weak_inlet.x, h=h)
                    T_points.append(state.T)
                    q = m_weak * (state.h - weak_inlet.h)
                    q_points.append(q)
                    x_points.append(state.x)
            else:
                vapor = weak_inlet

            # 2. Then cool to saturated liquid.
            for i, Qu in enumerate(np.linspace(vapor.Qu, 0, 10, endpoint=True)):
                state = amm.props2(P=weak_inlet.P, x=weak_inlet.x, Qu=Qu)
                T_points.append(state.T)
                q = m_weak * (state.h - weak_inlet.h)
                q_points.append(q)
                x_points.append(weak_inlet.x)
            self.q_pre = q

            T_points.pop()
            q_points.pop()
            x_points.pop()

        else:
            self.sat_inlet = weak_inlet
            self.m_sat = m_weak

        # Now we continue from the saturated liquid inlet towards the refrigerant.
        # TODO: need a robust way to choose xmax.
        # If xmax is less that sat_inlet.x, then the points will "go the wrong way".
        xmax = self.refrig_inlet.x * (0.8)
        for (i, x) in enumerate(np.linspace(self.sat_inlet.x, xmax)):
            try:
                q, t = self._x(x)
                x_points.append(x)
                q_points.append(q)
                T_points.append(t)
            except:
                print("[{}] x = {}: Unable to converge in NH3H2O".format(i, x))
                break

        if debug:
            print("Weak inlet: ", weak_inlet)
            print("Refrig inlet: ", refrig_inlet)
            print("Saturation point: ", self.sat_inlet)
            print(tabulate.tabulate(zip(x_points, T_points, q_points),
                                    ["x", "T", "Q"]))
        if debug:
            import matplotlib.pyplot as plt
            plt.figure()
            plt.title("Absorber stream, inlet Qu = {}".format(weak_inlet.Qu))
            plt.xlabel("Q")
            plt.ylabel("T")
            plt.plot(0, weak_inlet.T, 'o')
            plt.plot(q_points, T_points, '.')
            plt.show()

        self.q = scipy.interpolate.PchipInterpolator(T_points[::-1], q_points[::-1])
        self.T = scipy.interpolate.PchipInterpolator(q_points[::-1], T_points[::-1])

        if debug:
            q_range = np.linspace(min(q_points), max(q_points))
            T_vals = self.T(q_range)
            plt.plot(q_range, T_vals, '-')
            T_range = np.linspace(min(T_points), max(T_points))
            q_vals = self.q(T_range)
            plt.plot(q_vals, T_range, '--')
            plt.show()

    def _x(self, x_local):
        """ Determine the refrigerant absorbed and local solution mass flow rate.
        """
        x_weak = self.sat_inlet.x
        x_refrig = self.refrig_inlet.x
        x_rich = x_local
        m_weak = self.m_sat
        m_rich = m_weak * (x_weak - x_refrig) / (x_rich - x_refrig)
        m_refrig = m_weak * (x_weak - x_rich) / (x_rich - x_refrig)

        local_state = amm.props2(P=self.weak_inlet.P, Qu=0, x=x_local)
        # print(local_state)
        Q = self.q_pre + m_rich * local_state.h - m_weak * self.sat_inlet.h \
            - m_refrig * self.refrig_inlet.h
        # print(m_refrig,m_rich,Q,local_state.T)
        return Q, local_state.T


class AmmoniaGeneratorStream(object):
    """TODO: Account for non-saturated solution inlet.
        
    Inputs
    ------
    rich_inlet: ammonia_props.State
        The fluid state at the rich inlet.
    m_rich: number (kg/s)
        Inbound mass flow rate at the rich solution inlet.
    reflux_inlet: ammonia_props.State
        The fluid state at the liquid reflux inlet.
    m_reflux: number (kg/s)
        Inbound mass flow rate at the liquid reflux inlet.
    vapor_outlet: ammonia_props.State
        Fluid state at the vapor outlet.
    m_vapor: number (kg/s)
        Assumed outbound mass flow rate at the vapor outlet.
    """

    def __init__(self, rich_inlet, m_rich, reflux_inlet, m_reflux, vapor_outlet, m_vapor, debug=False):

        self.rich_inlet = rich_inlet
        self.m_rich = m_rich
        self.reflux_inlet = reflux_inlet
        self.m_reflux = m_reflux
        self.vapor_outlet = vapor_outlet
        self.m_vapor = m_vapor
        # Weak solution outlet is known...
        self.m_net = self.m_rich + self.m_reflux - self.m_vapor
        m_ammonia_net = self.rich_inlet.x * self.m_rich \
                        + self.reflux_inlet.x * self.m_reflux \
                        - self.vapor_outlet.x * self.m_vapor
        self.x_net = m_ammonia_net / self.m_net

        T_points = []
        q_points = []
        x_points = []
        self.q_pre = 0

        if self.rich_inlet.isSubcooled():
            # Need to heat it up.
            liquid, vapor = amm.equilibriumStates(self.rich_inlet.P, rich_inlet.x)
            # TODO: This is not graceful, but endpoint is used for the final
            # quantities. Eg., q_pre is used in self._x()
            for h in np.linspace(rich_inlet.h, liquid.h, 10, endpoint=True):
                state = amm.props2(P=rich_inlet.P, x=rich_inlet.x, h=h)
                q = m_rich * (h - rich_inlet.h)
                q_points.append(q)
                T_points.append(state.T)
                x_points.append(state.x)
            self.q_pre = q
            self.sat_inlet = liquid
            self.m_sat = m_rich
            self.m_desorb = m_vapor

            q_points.pop()
            T_points.pop()
            x_points.pop()

        elif self.rich_inlet.isSuperheated():
            # By design, this should not happen! Just use a rectifier.
            raise ValueError("Generator inlet was superheated.")
        else:
            # Already have some vapor, so we can separate it.
            # It will go into the rectifier. Assume:
            #   1. This flash vapor contributes to the specified outbound rate.
            liquid, vapor = amm.equilibriumStates(self.rich_inlet.P, rich_inlet.x)
            self.m_sat = m_rich * (1 - self.rich_inlet.Qu)
            self.sat_inlet = liquid
            self.m_desorb = m_vapor - m_rich * self.rich_inlet.Qu
            self.q_pre = 0
        
        # Generate a list of state points (and cum. heat flows)
        # by varying mass fraction from saturated inlet *down* to an arbitrary
        # low value (that should work to integrate heat exchange with any given
        # other stream).
        # Previously, I chose x_low_trace = 0.1
        # TODO: Need to assert self.sat_inlet.x > x_low_trace, or change limit.

        for (i, x) in enumerate(np.linspace(self.sat_inlet.x, 0.1)):
            # try:
            q, t = self._x(x)
            q_points.append(q)
            T_points.append(t)
            x_points.append(x)
            # except StandardError as e:
            #    print("[{}] x = {}: Unable to converge in NH3H2O".format(i,x))
            #    print(e)
            #    break

        if debug:
            print(tabulate.tabulate(zip(q_points, T_points, x_points),
                                    "q T x".split()))
            print(np.diff(T_points) < 0)

        self.q = scipy.interpolate.PchipInterpolator(T_points, q_points)
        self.T = scipy.interpolate.PchipInterpolator(q_points, T_points)

        if debug:
            import matplotlib.pyplot as plt
            plt.figure()
            plt.title("Generator stream, inlet Qu = {}".format(rich_inlet.Qu))
            plt.plot(self.q_pre, self.sat_inlet.T, 'o')
            plt.plot(q_points, T_points, '.')
            q_range = np.linspace(min(q_points), max(q_points))
            T_vals = self.T(q_range)
            plt.plot(q_range, T_vals, '-')
            T_range = np.linspace(min(T_points), max(T_points))
            q_vals = self.q(T_range)
            plt.plot(q_vals, T_range, '--')
            plt.show()
            print(tabulate.tabulate(zip(q_points, T_points, x_points),
                                    "q T x".split()))
            print(np.diff(T_points) < 0)

    def _x(self, z_local):
        # Input z_local is the ammonia mass fraction in the liquid phase.
        liquid, vapor = amm.equilibriumStates(self.rich_inlet.P, z_local)
        # Determine the amount of refrigerant vaporized and local states.
        x_net = self.x_net
        m_net = self.m_net
        x_liquid = z_local
        x_vapor = vapor.x

        m_liquid = m_net * (x_net - x_vapor) / (x_liquid - x_vapor)
        m_vapor = m_net * (x_net - x_liquid) / (x_liquid - x_vapor)

        Q = m_liquid * liquid.h - m_vapor * vapor.h \
            - self.m_rich * self.rich_inlet.h \
            - self.m_reflux * self.reflux_inlet.h \
            + self.m_vapor * self.vapor_outlet.h

        return Q, liquid.T


class AmmoniaRefluxStream(object):
    """A class to provide heat curves for the reflux coil.
    
    The model assumes internal vapor-liquid equilibrium and counterflow,
    and is therefore nearly identical to generator model, except that the
    generator has additional ports.
    
    For sign convention to work, we have to use the bottom of the rectifier
    (liquid outlet, vapor inlet) as the stream inlet where Q=0, since Q < 0
    as the process proceeds and T_top < T_bottom. Thus, since we would have to
    completely specify the process anyway, we might as well input the conditions
    at this point in the process. Driving the reaction past the point
    where liquid flow rate is zero might be a bad idea.
    
    Inputs
    ------
    vapor_inlet: ammonia_props.State
        The fluid state at the (medium ammonia) vapor inlet.
    m_inlet: number (kg/s)
        Inbound mass flow rate at the (medium ammonia) vapor inlet.
    liquid_outlet: ammonia_props.State
        The fluid state at the liquid reflux (low ammonia) outlet.
    m_reflux: number (kg/s)
        Outbound mass flow rate at the liquid reflux outlet.
    """

    def __init__(self, vapor_inlet, m_inlet, liquid_outlet, m_reflux, debug=False):
        self.vapor_inlet = vapor_inlet
        self.m_inlet = m_inlet
        self.liquid_outlet = liquid_outlet
        self.m_reflux = m_reflux

        # Determine net flow rate and net ammonia fraction
        self.m_net = self.m_inlet - self.m_reflux
        self.ammonia_net = m_inlet * vapor_inlet.x - m_reflux * liquid_outlet.x
        good = (self.m_net >= 0) and (self.ammonia_net >= 0)
        if not good:
            raise ValueError("In rectifier, net mass or ammonia flow is negative.")
        self.x_net = self.ammonia_net / self.m_net

        q_points, T_points = [], []
        #a_T_vapor, a_x_liquid, a_h_liquid, a_h_vapor = None, None, None, None
        list_vapor, list_liquid = None, None
        if debug:
            #a_T_vapor, a_x_liquid, a_h_liquid, a_h_vapor = [], [], [], []
            list_vapor, list_liquid = [], []

            print("=" * 20 + "Rectifier debug info" + "=" * 20)
            print("vapor_inlet (gen_vapor_outlet) = ", self.vapor_inlet)
            print("liquid_outlet (gen_reflux_inlet) = ", self.liquid_outlet)
            print("m_net, x_net = ", self.m_net, self.x_net)

        # This lookup has a bug: lookups are a non-monotonic function of x.
        #z_points = np.linspace(self.x_net, vapor_inlet.x)
        #for z in z_points:
        #    q, T = self._x(z, list_vapor, list_liquid)
        #    q_points.append(q)
        #    T_points.append(T)

        # Instead, let's parameterize by T.
        # Figure out the minimum temperature we may need to estimate.
        # Start by a lookup with inputs (P, x, Qu) to estimate T,
        # then switch and use inputs (P, T, Qu) with iterative solve to match x.
        T_guess = amm.props2(P=self.vapor_inlet.P, x=self.x_net, Qu=1).T
        x_error = lambda alpha: (self.x_net - amm.props2(
            P=self.vapor_inlet.P,
            T=T_guess + np.exp(alpha),
            Qu=1).x)
        if debug:
            print("T_guess (x = {}) = {}".format(self.x_net, T_guess))

        #opt = scipy.optimize.fsolve(x_error, -100.)[0]
        opt = -100.
        T_min = T_guess + np.exp(opt)
        T_max = vapor_inlet.T
        if debug:
            print("T_min, T_max = {}, {}".format(T_min, T_max))
        T_points = np.linspace(T_min, T_max)
        for t in T_points:
            q, _ = self._t(t, list_vapor, list_liquid)
            q_points.append(q)
        q_points = np.array(q_points)

        if debug:
            qdiff = np.diff(q_points)
            print(qdiff < 0)

            a_vapor = convert_state_list_to_array(list_vapor)
            a_liquid = convert_state_list_to_array(list_liquid)

            print(tabulate.tabulate(
                zip(a_vapor['T'], a_liquid['T'], a_vapor['x'], a_liquid['x'],
                    a_vapor['h'], a_liquid['h'], q_points),
                ["T vapor", "T liquid", "x vapor", "x liquid", "h vapor",
                 "h liquid", "q (kW)"]))

        # TODO: check for non-increasing points before interpolate.
        # Check for crossing over saturation and other causes...

        self.q = scipy.interpolate.PchipInterpolator(T_points, q_points)
        self.T = scipy.interpolate.PchipInterpolator(q_points, T_points)


    def _x(self, z_local, output_vapor=None, output_liquid=None):
        """Returns the amount of heat removed from the reflux stream between the
        cross section where vapor enters and the cross section given, subject
        to mass and species flow conservation and equilibrium.

        Input z_local is the ammonia mass fraction in the vapor phase.

        Warning: The underlying property function is not monotonic! Avoid this.
        """

        liquid, vapor = amm.equilibriumStates2(self.vapor_inlet.P, z_local)
        # print(liquid, vapor)
        # Determine the amount of refrigerant vaporized and local states.
        x_net = self.x_net
        m_net = self.m_net
        x_vapor = z_local
        x_liquid = liquid.x

        m_liquid = m_net * (x_vapor - x_net) / (x_liquid - x_vapor)
        m_vapor = m_net * (x_liquid - x_net) / (x_liquid - x_vapor)

        Q = m_vapor * vapor.h - m_liquid * liquid.h \
            + self.m_reflux * self.liquid_outlet.h \
            - self.m_inlet * self.vapor_inlet.h

        if output_vapor is not None:
            output_vapor.append(vapor)
        if output_liquid is not None:
            output_liquid.append(liquid)

        return Q, liquid.T

    def _t(self, t_local, output_vapor=None, output_liquid=None):
        """Returns the amount of heat removed from the reflux stream between the
        cross section where vapor enters and the cross section given, subject
        to mass and species flow conservation and equilibrium.

        Input t_local is the local temperature.
        """

        liquid, vapor = amm.equilibriumStates3(self.vapor_inlet.P, t_local)
        # Determine the amount of refrigerant vaporized and local states.
        x_net = self.x_net
        m_net = self.m_net
        x_vapor = vapor.x
        x_liquid = liquid.x

        m_liquid = m_net * (x_vapor - x_net) / (x_liquid - x_vapor)
        m_vapor = m_net * (x_liquid - x_net) / (x_liquid - x_vapor)

        Q = m_vapor * vapor.h - m_liquid * liquid.h \
            + self.m_reflux * self.liquid_outlet.h \
            - self.m_inlet * self.vapor_inlet.h

        if output_vapor is not None:
            output_vapor.append(vapor)
        if output_liquid is not None:
            output_liquid.append(liquid)

        return Q, liquid.T


class AmmoniaChiller(object):
    def __init__(self):
        self.points = """rich abs outlet
rich pump outlet
rich shx outlet
rich gen sat liquid
weak gen outlet
weak shx outlet
weak exp outlet
gen vapor outlet
gen reflux inlet
refrig rect outlet
refrig cond outlet
refrig cehx liquid outlet
refrig exp outlet
refrig evap outlet
refrig cehx sat vapor
refrig cehx vapor outlet
rectifier_liquid
gen_vapor_formation
abs_vapor_final""".replace(" ", "_").split('\n')
        vars = """Q_abs,kW
Q_gen,kW
Q_cond,kW
Q_evap,kW
Q_reflux,kW
Q_shx,kW
Q_cehx,kW
W_pump,kW
COP,kW/kW
ZeroCheck,kW
ZeroCheckSHX,kW
ZeroCheckCEHX,kW
Q_gen_rect_combo,kW
Q_refrig_side,kW
ZeroCheckRect,kg/s
m_rich,kg/s
m_weak,kg/s
m_gen_vapor,kg/s
m_gen_reflux,kg/s
m_refrig,kg/s
check_rectifier_delta_T,K
is_degenerate,bool""".split()
        self.vars = [var.split(',')[0] for var in vars]
        self.units = [var.split(',')[1] for var in vars]

        self.rich_abs_outlet = amm.props2(T=400, P=10, x=0.5)
        self.rich_pump_outlet = amm.props2(T=400, P=10, x=0.5)
        self.rich_shx_outlet = amm.props2(T=400, P=10, x=0.5)
        self.rich_gen_sat_liquid = amm.props2(T=400, P=10, x=0.5)
        self.weak_gen_outlet = amm.props2(T=400, P=10, x=0.5)
        self.weak_shx_outlet = amm.props2(T=400, P=10, x=0.5)
        self.weak_exp_outlet = amm.props2(T=400, P=10, x=0.5)
        self.gen_vapor_outlet = amm.props2(T=400, P=10, x=0.5)
        self.gen_reflux_inlet = amm.props2(T=400, P=10, x=0.5)
        self.refrig_rect_outlet = amm.props2(T=400, P=10, x=0.5)
        self.refrig_cond_outlet = amm.props2(T=400, P=10, x=0.5)
        self.refrig_cehx_liquid_outlet = amm.props2(T=400, P=10, x=0.5)
        self.refrig_exp_outlet = amm.props2(T=400, P=10, x=0.5)
        self.refrig_evap_outlet = amm.props2(T=400, P=10, x=0.5)
        self.refrig_cehx_sat_vapor = amm.props2(T=400, P=10, x=0.5)
        self.refrig_cehx_vapor_outlet = amm.props2(T=400, P=10, x=0.5)
        self.rectifier_liquid = amm.props2(T=400, P=10, x=0.5)
        self.gen_vapor_formation = amm.props2(T=400, P=10, x=0.5)
        self.abs_vapor_final = amm.props2(T=400, P=10, x=0.5)

        self.Q_abs = 0
        self.Q_gen = 0
        self.Q_cond = 0
        self.Q_evap = 0
        self.Q_reflux = 0
        self.Q_shx = 0
        self.Q_cehx = 0
        self.COP = 0
        self.ZeroCheck = 0
        self.W_pump = 0
        self.m_rich = 1
        self.m_weak = 0
        self.m_gen_vapor = 0
        self.m_gen_reflux = 0
        self.m_refrig = 0
        self.ZeroCheckSHX = 0
        self.ZeroCheckCEHX = 0
        self.Q_gen_rect_combo = 0
        self.Q_refrig_side = 0
        self.ZeroCheckRect = 0
        self.check_rectifier_delta_T = 0
        self.is_degenerate = False

    def getStateIterator(self):
        return stateIterator(self)

    def stateTable(self):
        ii = range(len(self.points))
        states = np.zeros_like(ii, dtype=StateType)
        for i, s2 in zip(ii, stateIterator(self)):
            states[i] = s2
        return states

    def getStateTable(self):
        return CStateTable(self.stateTable(), self.points)

    def getVariablesTable(self):
        vals = [self.__getattribute__(var) for var in self.vars]
        return CVariablesTable(self.vars, vals, self.units)

    def __repr__(self):
        return "{}\n{}".format(
            self.getStateTable().__repr__(),
            self.getVariablesTable().__repr__()
        )

    def _repr_html_(self):
        template = """<h3>State points</h3>
        {}
        <br/>
        <h3>Performance variables</h3>
        {}
        """
        return template.format(
            self.getStateTable()._repr_html_(),
            self.getVariablesTable()._repr_html_()
        )

    def update(self,
               x_refrig=0.999869,
               T_evap=5.3 + 273.15,
               T_cond=38.7 + 273.15,
               Qu_evap=0.998,
               eff_CEHX=0.95,
               T_abs_outlet=37 + 273.15,
               T_gen_outlet=101 + 273.15,
               m_rich=0.40,
               eta_pump=0.8,
               eff_SHX=0.7,
               T_rect=40.5 + 273.15):
        """Solves the cycle.
        
        Refrigerant cycle step 1:
        1. Input condenser and evaporator outlet temperatures, mass fraction.
        2. Compute condenser and evaporator pressures.
        Refrigerant CEHX:
        1. Input inlet states, effectiveness.
        2. Compute outlet states.
        Refrigerant expansion valve:
        1. Input inlet state
        2. Compute outlet state
        Absorber step 1:
        1. Input absorber pressure and outlet temperature.
        2. Compute outlet mass fraction.
        Generator step 1:
        1. Input generator pressure and outlet temperature.
        2. Compute outlet mass fraction.
        Refrigerant cycle step 2:
        1. Compute mass flow rates.
        Pump:
        1. Input inlet state and outlet pressure.
        2. Compute outlet state.
        SHX:
        1. Input inlet states and effectiveness.
        2. Output outlet states.
        Solution expansion valve:
        1. Input inlet state.
        2. Compute outlet state.
        Absorber step 2:
        1. Input inlet temperature.
        2. Determine heat flow and saturation temperature.
        Generator step 2:
        1. Input inlet temperature.
        2. Determine heat flow and saturation temperature.

        TODO: figure out how T_rect is a relevant input.
        """

        self.is_degenerate = False

        # Refrigerant cycle step 1:
        self.refrig_evap_outlet, self.refrig_cond_outlet \
            = self.updateRefrig(x_refrig, T_evap, T_cond, Qu_evap)

        # Refrigerant CEHX:
        self.refrig_cehx_liquid_outlet, \
        self.refrig_cehx_sat_vapor, \
        self.refrig_cehx_vapor_outlet \
            = self.updateCEHX(self.refrig_cond_outlet,
                              self.refrig_evap_outlet,
                              eff_CEHX)

        # Refrigerant expansion valve:
        self.refrig_exp_outlet \
            = self.updateExpander(self.refrig_cehx_liquid_outlet,
                                  self.refrig_evap_outlet.P)

        # Absorber step 1:
        self.rich_abs_outlet \
            = self.updateAbsorber1(self.refrig_cehx_vapor_outlet.P,
                                   T_abs_outlet)

        # Generator step 1:
        self.weak_gen_outlet \
            = self.updateGenerator1(self.refrig_cond_outlet.P,
                                    T_gen_outlet)

        # Refrigerant cycle step 2:
        self.m_rich = m_rich
        x_weak, x_rich = self.weak_gen_outlet.x, self.rich_abs_outlet.x
        self.m_refrig, self.m_weak = self.updateFlowRates(
            m_rich, x_rich, x_weak, x_refrig)

        # Pump:
        self.rich_pump_outlet = self.updatePump(self.rich_abs_outlet,
                                                self.refrig_cond_outlet.P,
                                                eta_pump)

        # SHX:
        self.rich_shx_outlet, self.weak_shx_outlet = \
            self.updateSHX(self.rich_pump_outlet, self.weak_gen_outlet,
                           eff_SHX, self.m_rich, self.m_weak)

        # Solution expansion valve:
        self.weak_exp_outlet = \
            self.updateExpander(self.weak_shx_outlet, self.refrig_evap_outlet.P)

        # Update some saturated states
        self.rich_gen_sat_liquid = amm.props2(x=self.rich_shx_outlet.x,
                                              P=self.rich_shx_outlet.P,
                                              Qu=0)
        self.gen_vapor_outlet = amm.props2(P=self.rich_gen_sat_liquid.P,
                                           T=self.rich_gen_sat_liquid.T,
                                           Qu=1)
        self.gen_reflux_inlet = amm.props2(P=self.rich_gen_sat_liquid.P,
                                           T=self.rich_gen_sat_liquid.T,
                                           Qu=0)
        self.refrig_rect_outlet = amm.props2(P=self.refrig_cond_outlet.P,
                                             x=x_refrig,
                                             T=T_rect)
        # Should be positive
        self.check_rectifier_delta_T = self.gen_vapor_outlet.T - T_rect
        self.is_degenerate = self.is_degenerate or self.check_rectifier_delta_T < 0

        x_vapor, x_liquid = self.gen_vapor_outlet.x, self.gen_reflux_inlet.x
        self.m_gen_vapor, self.m_gen_reflux = \
            self.updateRectifier(self.m_refrig, x_refrig, x_vapor, x_liquid)
        self.rectifier_liquid = amm.props2(T=T_rect,
                                           P=self.refrig_rect_outlet.P,
                                           Qu=0)
        self.gen_vapor_formation = amm.props2(T=T_gen_outlet,
                                              P=self.weak_gen_outlet.P,
                                              Qu=1)
        self.abs_vapor_final = amm.props2(T=self.weak_exp_outlet.T,
                                          P=self.weak_exp_outlet.P,
                                          Qu=1)

        self.Q_cond = self.m_refrig * self.refrig_cond_outlet.h \
                      - self.m_refrig * self.refrig_rect_outlet.h
        self.Q_evap = self.m_refrig * self.refrig_evap_outlet.h \
                      - self.m_refrig * self.refrig_exp_outlet.h

        self.updateAbsorber2()
        self.updateGenerator2()

        self.Q_reflux = self.m_refrig * self.refrig_rect_outlet.h \
                        + self.m_gen_reflux * self.gen_reflux_inlet.h \
                        - self.m_gen_vapor * self.gen_vapor_outlet.h
        self.W_pump = self.m_rich * self.rich_pump_outlet.h \
                      - self.m_rich * self.rich_abs_outlet.h
        self.ZeroCheck = self.Q_abs + self.Q_cond + self.Q_gen \
                         + self.Q_reflux + self.Q_evap + self.W_pump
        self.COP = self.Q_evap / self.Q_gen

        # More debugging
        self.Q_SHXhot = self.m_weak * self.weak_shx_outlet.h \
                        - self.m_weak * self.weak_gen_outlet.h
        self.Q_shx = self.m_rich * self.rich_shx_outlet.h \
                     - self.m_rich * self.rich_pump_outlet.h
        self.ZeroCheckSHX = self.Q_shx + self.Q_SHXhot
        self.Q_CEHXhot = self.m_refrig * self.refrig_cehx_liquid_outlet.h \
                         - self.m_refrig * self.refrig_cond_outlet.h
        self.Q_cehx = self.m_refrig * self.refrig_cehx_vapor_outlet.h \
                      - self.m_refrig * self.refrig_evap_outlet.h
        self.ZeroCheckCEHX = self.Q_CEHXhot + self.Q_cehx

        self.Q_gen_rect_combo = self.m_refrig * self.refrig_rect_outlet.h \
                                + self.m_weak * self.weak_gen_outlet.h \
                                - self.m_rich * self.rich_shx_outlet.h
        self.Q_refrig_side = self.m_refrig * self.refrig_rect_outlet.h \
                             - self.m_refrig * self.refrig_cehx_vapor_outlet.h
        self.ZeroCheckRect = self.m_refrig * x_refrig \
                             + self.m_gen_reflux * x_liquid \
                             - self.m_gen_vapor * x_vapor

    def updateRefrig(self, x_refrig, T_evap, T_cond, Qu_evap):
        evap_outlet = amm.props2(T=T_evap,
                                 x=x_refrig,
                                 Qu=Qu_evap)
        cond_outlet = amm.props2(T=T_cond,
                                 x=x_refrig,
                                 Qu=0)
        return evap_outlet, cond_outlet

    def updateCEHX(self, liquid_inlet, vapor_inlet, effectiveness):
        liquid_max = amm.props2(T=vapor_inlet.T,
                                x=liquid_inlet.x,
                                P=liquid_inlet.P)
        vapor_max = amm.props2(T=liquid_inlet.T,
                               x=vapor_inlet.x,
                               P=vapor_inlet.P)
        deltaH_liquid_max = liquid_inlet.h - liquid_max.h
        deltaH_vapor_max = vapor_max.h - vapor_inlet.h
        deltaH_max = min(deltaH_liquid_max, deltaH_vapor_max)
        deltaH = effectiveness * deltaH_max
        liquid_outlet = amm.props2(x=liquid_inlet.x,
                                   P=liquid_inlet.P,
                                   h=liquid_inlet.h - deltaH)
        vapor_outlet = amm.props2(x=vapor_inlet.x,
                                  P=vapor_inlet.P,
                                  h=vapor_inlet.h + deltaH)
        sat = amm.props2(x=vapor_inlet.x,
                         P=vapor_inlet.P,
                         Qu=1)
        return liquid_outlet, sat, vapor_outlet

    def updateExpander(self, inlet, P_outlet):
        outlet = amm.props2(x=inlet.x, h=inlet.h, P=P_outlet)
        return outlet

    def updateAbsorber1(self, P_inlet, T_outlet):
        outlet = amm.props2(P=P_inlet, T=T_outlet, Qu=0)
        return outlet

    def updateGenerator1(self, P_inlet, T_outlet):
        outlet = amm.props2(P=P_inlet, T=T_outlet, Qu=0)
        return outlet

    def updateFlowRates(self, m_rich, x_rich, x_weak, x_refrig):
        # m_rich = m_weak + m_refrig
        # m_rich * x_rich = m_refrig * x_refrig + m_weak * x_weak
        m_weak = m_rich * (x_rich - x_refrig) / (x_weak - x_refrig)
        m_refrig = m_rich * (x_rich - x_weak) / (x_refrig - x_weak)
        return m_refrig, m_weak

    def updatePump(self, inlet, P_outlet, eta_pump):
        """Pump:
        1. Input inlet state and outlet pressure.
        2. Compute outlet state."""
        outlet_ideal = amm.props2(x=inlet.x, P=P_outlet, s=inlet.s)
        deltaH_ideal = outlet_ideal.h - inlet.h
        deltaH = deltaH_ideal / eta_pump
        h_out = inlet.h + deltaH
        outlet = amm.props2(x=inlet.x, P=P_outlet, h=h_out)
        return outlet

    def updateSHX(self, cold_inlet, hot_inlet, effectiveness, m_cold, m_hot):
        """SHX:
        1. Input inlet states and effectiveness.
        2. Output outlet states."""
        cold_outlet = amm.props2(T=hot_inlet.T,
                                 x=cold_inlet.x,
                                 P=cold_inlet.P)
        hot_outlet = amm.props2(T=cold_inlet.T,
                                x=hot_inlet.x,
                                P=hot_inlet.P)
        deltaH_cold_max = cold_outlet.h - cold_inlet.h
        deltaH_hot_max = hot_inlet.h - hot_outlet.h
        Q_max = min(m_cold * deltaH_cold_max, m_hot * deltaH_hot_max)
        Q = effectiveness * Q_max
        cold_outlet = amm.props2(x=cold_inlet.x,
                                 P=cold_inlet.P,
                                 h=cold_inlet.h + Q / m_cold)
        hot_outlet = amm.props2(x=hot_inlet.x,
                                P=hot_inlet.P,
                                h=hot_inlet.h - Q / m_hot)
        # Cold/Hot refer to the streams
        return cold_outlet, hot_outlet

    def updateAbsorber2(self):
        """Absorber step 2:
        1. Input inlet temperature.
        2. Determine heat flow and saturation temperature.
        
        Absorber outlet streams:
            Ammonia-rich solution to generator via pump and SHX
        Absorber inlet streams: 
            Ammonia-weak solution from generator via SHX & expander
            (Mostly vapor) two-state refrigerant from CEHX
        """
        self.Q_abs = self.m_rich * self.rich_abs_outlet.h \
                     - self.m_weak * self.weak_exp_outlet.h \
                     - self.m_refrig * self.refrig_cehx_vapor_outlet.h

    def updateGenerator2(self):
        """Generator step 2:
        1. Input inlet temperature.
        2. Determine heat flow and saturation temperature.
        
        Generator outlet streams:
            Ammonia-weak solution to absorber via SHX & expander
            (Almost Ammonia) Generator vapor outlet
        Generator inlet streams:        
            Ammonia-rich solution from absorber via pump & SHX
            Ammonia-rich solution return from rectifier
        """
        self.Q_gen = self.m_weak * self.weak_gen_outlet.h \
                     + self.m_gen_vapor * self.gen_vapor_outlet.h \
                     - self.m_rich * self.rich_shx_outlet.h \
                     - self.m_gen_reflux * self.gen_reflux_inlet.h

    def updateRectifier(self, m_refrig, x_refrig, x_vapor, x_liquid):
        """Return the mass flow rates of vapor leaving generator and liquid
        returning to generator, respectively, at the rectifier interface."""
        # m_refrig = m_vapor - m_liquid
        # m_refrig * x_refrig = m_vapor * x_vapor - m_liquid * x_liquid
        m_vapor = m_refrig * (x_refrig - x_liquid) / (x_vapor - x_liquid)
        m_liquid = m_refrig * (x_refrig - x_vapor) / (x_vapor - x_liquid)
        return m_vapor, m_liquid

    def getPaths(self):
        return [(0, 1),
                (1, 2),
                (2, 3),
                (3, 7),
                (7, 8),
                (8, 3),
                (7, 9),
                (9, 16),
                (16, 8),
                (3, 4),
                (4, 17),
                (17, 7),
                (4, 5),
                (5, 6),
                (6, 0),
                (9, 10),
                (10, 11),
                (11, 12),
                (12, 13),
                (13, 14),
                (14, 15),
                (15, 0),
                (15, 18),
                (18, 6)]

    def getAbsorberStream(self):
        return AmmoniaAbsorberStream(self.weak_exp_outlet,
                                     self.refrig_cehx_vapor_outlet,
                                     self.m_weak)

    def getGeneratorStream(self, retry=True):
        args = (self.rich_shx_outlet, self.m_rich,
                self.gen_reflux_inlet, self.m_gen_reflux,
                self.gen_vapor_outlet, self.m_gen_vapor)
        try:
            return AmmoniaGeneratorStream(*args)
        except Exception as e:
            if retry:
                return AmmoniaGeneratorStream(*args, debug=True)
            else:
                raise e

    def getEvaporatorStream(self):
        return HRHX_integral_model.aquaStream(self.refrig_exp_outlet,
                                              self.m_refrig)

    def getCondenserStream(self):
        return HRHX_integral_model.aquaStream(self.refrig_rect_outlet,
                                              self.m_refrig)

    def getRectifierStream(self,retry=True):
        try:
            return AmmoniaRefluxStream(self.gen_vapor_outlet, self.m_gen_vapor,
                                         self.gen_reflux_inlet, self.m_gen_reflux)
        except Exception as e:
            if retry:
                return AmmoniaRefluxStream(
                    self.gen_vapor_outlet, self.m_gen_vapor,
                    self.gen_reflux_inlet, self.m_gen_reflux,
                    debug=True
                )
            else:
                raise e

    def getSHX(self, **kwargs):
        """Constructs and returns a model for the internal heat
        exchanger between the weak and rich solution streams.
        """
        hot = HRHX_integral_model.aquaStream(self.weak_gen_outlet, self.m_weak)
        cold = HRHX_integral_model.aquaStream(self.rich_pump_outlet, self.m_rich)
        return HRHX_integral_model.counterflow_integrator(cold, hot, **kwargs)

    def getCEHX(self, **kwargs):
        """Constructs and returns a model for the internal heat
        exchanger between the condenser and evaporator.
        """
        hot = HRHX_integral_model.aquaStream(self.refrig_cond_outlet, self.m_refrig)
        cold = HRHX_integral_model.aquaStream(self.refrig_evap_outlet, self.m_refrig)
        return HRHX_integral_model.counterflow_integrator(cold, hot, **kwargs)

    def display(a):
        import matplotlib.pyplot as plt
        plt.figure()
        plt.title("Internal streams")
        absorber = a.getAbsorberStream()
        plotStream(absorber,
                   (1.05 * a.Q_abs, 0),
                   (a.weak_exp_outlet.T, a.rich_abs_outlet.T),
                   {"label":"absorber"})
        gen = a.getGeneratorStream()
        plotStream(gen,
                   (0, 1.05 * a.Q_gen),
                   (a.rich_shx_outlet.T, a.weak_gen_outlet.T),
                   {"label":"generator"})
        evap = a.getEvaporatorStream()
        plotStream(evap,
                   (0, 1.05 * a.Q_evap),
                   (a.refrig_exp_outlet.T, a.refrig_evap_outlet.T),
                   {"label":"evaporator"})
        cond = a.getCondenserStream()
        plotStream(cond,
                   (1.05 * a.Q_cond, 0),
                   (a.refrig_rect_outlet.T, a.refrig_cond_outlet.T),
                   {"label":"condenser"})
        rect = a.getRectifierStream()
        plotStream(rect,
                   (1.05 * a.Q_reflux, 0),
                   (a.refrig_rect_outlet.T, a.gen_vapor_outlet.T),
                   {"label":"rectifier"})
        plt.legend(loc='best')
        plt.xlabel("Heat (kW)")
        plt.ylabel("Temperature (K)")
        plt.ylim([273,400])

        # Internal streams
        shx = a.getSHX(initQmax=True) #, qmaxopts={"brute": True})
        HRHX_integral_model.plotFlow(shx, Qactual=a.Q_shx)
        plt.title("Solution heat exchanger (SHX)")
        plt.xlabel("Heat (kW)")
        plt.ylabel("Temperature (K)")
        cehx = a.getCEHX(initQmax=True) #, qmaxopts={"brute": True})
        HRHX_integral_model.plotFlow(cehx, Qactual=a.Q_cehx)
        plt.title("Condenser-Evaporator heat exchanger (CEHX)")
        plt.xlabel("Heat (kW)")
        plt.ylabel("Temperature (K)")


def plotStream(stream, qrange, Trange, plotopts={}):
    import matplotlib.pyplot as plt
    # First let's plot T(q) as line
    q1 = np.linspace(*qrange)
    T1 = stream.T(q1)
    l = plt.plot(q1, T1, '-', **plotopts)[0]
    # Now we print q(T) as dots
    T2 = np.linspace(*Trange)
    q2 = stream.q(T2)
    plt.plot(q2, T2, '.', color=l.get_c(), **plotopts)


def main():
    a = AmmoniaChiller()
    print(a)
    a.update()
    print(a)
    a.display()
    return a


if __name__ == "__main__":
    a = main()
