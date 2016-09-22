# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 08:21:01 2016

@author: nfette

Adsorption cycle simulation, fast optimization, and parametric study.
"""

from numpy import linspace, zeros, meshgrid, nan, logical_and, unravel_index
from scipy.interpolate import PchipInterpolator
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, plot, show, xlabel, ylabel, contour, clabel, title, subplot, xlim, draw
from mpl_toolkits.mplot3d import Axes3D
from adsorption import AdsorptionChiller, AdsorptionChillerControl, AdsorptionChillerSpec, WaterTQ2H

class convergeMe:
    """Adsorption cycle simulation and fast optimization

This routine effectively performs optimization of the dwell times in the
adsorption and desorption steps to maximize cooling capacity. Each cycle in an
adsorber includes:
    Desorption, decompression, adsorption, compression.

Compared to Gupta's original solution method, this solution:

* Adds dwell times for compression and decompression. These are significant
  fractions of total cycle time and should not be neglected. Computing
  them does not require integration.
* Integrates the ODE to find time given temperature, instead of finding
  temperature given a time interval. Temperature is a more universal input
  than cycle time. This applies only for this equilibrium model.
  
The optimization is a rough, first pass, intended to visualize the solution
space and the relations:
    
* between independent sets of internal control parameters, namely {adsorber bed
  temperature limits} vs. {dwell times}, and
* between internal control parameters {dwell times} and outputs {Q, COP}.

Numerically, the optimization uses a brute force comparison over the search
space which is a grid of points (q_low, q_high). The dwell times for adsorption
and desorption steps are interpolated from the discrete curves of adsorbed
ratio vs time, precomputed for those processes. Then, the compression and
decompression dwell times are computed at each point, as well as the resulting
time-averaged heat flow rates and COP.
"""

    def __init__(self,ch):
        t_cond = ch.ctrl.t_cond
        t_evap = ch.ctrl.t_evap
        t_cool = ch.ctrl.t_cool
        t_exhaust = ch.ctrl.t_exhaust
        
        # The extreme ranges depend on the heat and cold sources.
        self.q_min = ch.f.Q(t_cond, t_exhaust);
        self.q_max = ch.f.Q(t_evap, t_cool);
        q_range = linspace(self.q_min+0.0001,self.q_max-0.01);
        # We should start desorption at q_max.
        T_d0 = ch.f.T(t_cond, self.q_max-0.001);
        T_d1 = ch.f.T(t_cond, self.q_min+0.0001);
        
        T_a0 = ch.f.T(t_evap, self.q_min+0.0001);
        T_a1 = ch.f.T(t_evap, self.q_max-0.001);
        
        # Collect data for (all) desorption processes
        T_d = linspace(T_d0,T_d1);
        [self.t_d, self.q_d] = ch.desorptionFlip(T_d);
        
        # Collect data for (all) adsorption processes
        T_a = linspace(T_a0, T_a1);
        [self.t_a, self.q_a] = ch.adsorptionFlip(T_a);
        
        # Collect data for compression and decompression
        deltat_cmp = ch.compress(q_range)
        deltat_dec = ch.decompress(q_range)
        
        # Spline it.
        self.ppd = PchipInterpolator(self.q_d[::-1], self.t_d[::-1])
        self.ppa = PchipInterpolator(self.q_a, self.t_a)
        self.pp_comp = PchipInterpolator(q_range,deltat_cmp)
        self.pp_decomp = PchipInterpolator(q_range,deltat_dec)
        self.pp_Td = PchipInterpolator(self.q_d[::-1],T_d[::-1])
        self.pp_Ta = PchipInterpolator(self.q_a,T_a)
        
        # Stuff to display for user
        td_splined = self.ppd(q_range)
        ta_splined = self.ppa(q_range)
        
        figure(1)
        plt.cla()
        plot(self.t_d,self.q_d,'ro',
             self.t_a,self.q_a,'bo',
             td_splined,q_range,'r',
             ta_splined,q_range,'b')
        xlabel('Integration time, $t$ (s)')
        ylabel('Adsorbed mass ratio, $q$ (kg/kg)')
        title('$T_{{\\rm exhaust}}$ = {:g} K'.format(t_exhaust))
        
    def parametric(self,ch):
        """Now do a parametric study. Choose q_low and q_high, then determine
        the time needed in each step, and the other outputs of the system."""
        
        if False:
            Ni = 100
            Nj = 100
            self.q_lows = linspace(self.q_min + 0.0001, self.q_max - 0.001, Ni)
            self.q_highs = linspace(self.q_min + 0.0001, self.q_max - 0.001, Nj)
        else:
            self.q_lows = self.q_d[::1]
            self.q_highs = self.q_a[::1]
            Ni = len(self.q_lows)
            Nj = len(self.q_highs)
        self.t_desorb = zeros((Ni, Nj))
        self.t_decomp = zeros((Ni, Nj))
        self.t_adsorb = zeros((Ni, Nj))
        self.t_compress = zeros((Ni, Nj))
        
        self.t_cycle = zeros((Ni, Nj))
        self.q_ref = zeros((Ni, Nj))
        self.q_in = zeros((Ni, Nj))
        self.cop = zeros((Ni, Nj))
        mask1 = zeros((Ni, Nj))
        masknan = zeros((Ni, Nj))
        masknan.fill(nan)
        [X,Y] = meshgrid(self.q_lows,self.q_highs,indexing='ij')
        for i, q_low in enumerate(self.q_lows):
            self.t_desorb[i,:] += self.ppd(q_low)
            self.t_decomp[i,:] = self.pp_decomp(q_low)
            self.t_adsorb[i,:] += -self.ppa(q_low)
        for j, q_high in enumerate(self.q_highs):
            self.t_desorb[:,j] += - self.ppd(q_high)
            self.t_adsorb[:,j] += self.ppa(q_high)
            self.t_compress[:,j] = self.pp_comp(q_high)
        
        self.t_cycle = self.t_desorb + self.t_compress + self.t_adsorb + self.t_decomp
        #self.q_ref[i,j], self.q_in[i,j], self.cop[i,j], _ = ch.afterSolve(q_low,q_high,self.t_cycle[i,j])
        T_d1 = self.pp_Td(Y)
        T_a1 = self.pp_Ta(X)
        m_ref = 2 * (Y - X) * ch.spec.m2 / self.t_cycle
        DeltaH = WaterTQ2H(ch.ctrl.t_evap,1) - WaterTQ2H(ch.ctrl.t_cond,0)
        self.q_ref = m_ref * DeltaH
        self.q_in = (m_ref * ch.spec.hads) \
            + 2 * (ch.spec.m2 * (ch.spec.c2 + Y * ch.spec.cw) + ch.spec.m1 * ch.spec.c1) * (T_d1 - T_a1) / (self.t_cycle)
        self.cop = self.q_ref / self.q_in
        
        # Do not use this equation, it is wrong:
        #t_cycle = t_desorb + t_adsorb;
        # Desorb step fraction of the total time
        t_fractiond = self.t_desorb / self.t_cycle;
        mask1[logical_and(self.t_desorb >= 0, self.t_adsorb >= 0)] = 1;
        masknan[logical_and(self.t_desorb >= 0, self.t_adsorb >= 0)] = 1;
        Iflat = (self.q_ref * mask1).argmax()
        I1,I2 = unravel_index(Iflat,self.q_ref.shape)
        t_d_opt = self.t_desorb[I1,I2];
        t_e_opt = self.t_decomp[I1,I2];
        t_a_opt = self.t_adsorb[I1,I2];
        t_c_opt = self.t_compress[I1,I2];
        Q_max = self.q_ref[I1,I2];
        COP_max = self.cop[I1,I2];
        
        figure(2)
        plt.cla()
        #CS = contour(X,Y,t_cycle)
        #clabel(CS, inline=1, fontsize=10)
        plt.pcolormesh(X,Y,self.t_cycle)
        xlabel('Low mass ratio, $q_{low}$ (kg/kg)')
        ylabel('High mass ratio, $q_{high}$ (kg/kg)')
        title('Cycle time, $t_{cycle}$ (s)')
        
        fig=figure(3)
        ax = fig.gca()
        #CS = contour(X,Y,t_fractiond)
        #clabel(CS, inline=1, fontsize=10)
        #plt.pcolormesh(X,Y,masknan*t_fractiond,vmin=0,vmax=1)
        msh=plt.pcolormesh(X,Y,mask1*self.q_ref)
        xlabel('Low mass ratio, $q_{low}$ (kg/kg)')
        ylabel('High mass ratio, $q_{high}$ (kg/kg)')
        title('Desorption step fraction, $t_{des}/t_{cycle}$')
        
        fig = figure(4)
        plt.clf()
        ax = fig.add_subplot(111, projection='3d')
        ax.grid(True)
        ax.plot_wireframe(self.t_desorb*masknan,self.t_adsorb*masknan,X,color='b')
        ax.plot_wireframe(self.t_desorb*masknan,self.t_adsorb*masknan,Y,color='r')
        #plot3([t_d_opt t_d_opt], [t_a_opt t_a_opt], [X(I1,I2) Y(I1,I2)], 'ko-')
        #set(gca,'xlim',[0, max(t_desorb(:))])
        #set(gca,'ylim',[0, max(t_adsorb(:))])
        xlabel('Desorb step time, $\Delta t$ (s)')
        ylabel('Adsorb step time, $\Delta t$ (s)')
        ax.set_zlabel('Mass ratio, $q$ (kg/kg)')
        
        fig = figure(5)
        plt.clf()
        ax = fig.add_subplot(111, projection='3d')
        ax.grid(True)
        #plot3(t_desorb, t_adsorb, q_in, 'b-')
        #plot3(t_desorb', t_adsorb', q_in', 'r-')
        ax.plot_wireframe(self.t_desorb*masknan, self.t_adsorb*masknan, self.q_ref*masknan, color='b')
        ax.plot([t_d_opt, t_d_opt], [t_a_opt, t_a_opt], [0, Q_max], 'ko-')
        xlabel('Desorb step time, $\Delta t$ (s)')
        ylabel('Adsorb step time, $\Delta t$ (s)')
        ax.set_zlabel('Cooling capacity, $Q_{cool}$ (kW)')
            
        return [t_d_opt, t_e_opt, t_a_opt, t_c_opt, Q_max, COP_max]

class parstudy1():
    def __init__(self):
        #global t_exhaust T_exhaust_range
        spec = AdsorptionChillerSpec()
        ctrl = AdsorptionChillerControl()
        self.chiller = AdsorptionChiller(spec, ctrl)
        
        Ni = 5;
        #T_exhaust_range = linspace(313,393,Ni);
        self.T_exhaust_range = linspace(313,393,Ni)
        self.tt_d = zeros(Ni)
        self.tt_e = zeros(Ni)
        self.tt_a = zeros(Ni)
        self.tt_c = zeros(Ni)
        self.QQ = zeros(Ni)
        self.CC = zeros(Ni)
        
        figure(6)
        plt.clf()
        self.ax1=subplot(2,1,1)
        self.h = plot(self.T_exhaust_range-273,self.QQ,'ko')[0]
        xlim([40,120])
        plt.ylim([0,100])
        ylabel('$Q_{\\rm max}$ (kW)')
        self.ax2=subplot(2,1,2)
        self.h2 = plot(self.T_exhaust_range-273,self.CC,'ko')[0]
        xlim([40,120])
        plt.ylim([0,1])
        ylabel('COP')
        xlabel('Exhaust temperature, $T$ ($^\circ$C)')
        
        self.it = iter(enumerate(self.T_exhaust_range))
    
    def __iter__(self):
        return self
        
    def __next__(self):
        i,t_exhaust = self.it.__next__()
        
        self.chiller.ctrl.t_exhaust = t_exhaust
        print('{} K, '.format(t_exhaust))
        self.something = convergeMe(self.chiller)
        t_d_opt, t_e_opt, t_a_opt, t_c_opt, Q_max, COP_max \
            = self.something.parametric(self.chiller)
        self.tt_d[i] = t_d_opt
        self.tt_e[i] = t_e_opt
        self.tt_a[i] = t_a_opt
        self.tt_c[i] = t_c_opt
        self.QQ[i] = Q_max
        self.CC[i] = COP_max
        print('{:g} s, {:g} s, {:g} s, {:g} s, {:g} kW\n'.format(t_d_opt, t_e_opt, t_a_opt, t_c_opt, Q_max))
        self.h.set_ydata(self.QQ)
        self.h2.set_ydata(self.CC)
        for f in range(1,7):
            fig = figure(f)
            fig.canvas.draw()
            fig.canvas.flush_events()

if __name__ == "__main__":
    #input('Hit enter to proceed')
    ps = parstudy1()
    for p in ps:
        pass