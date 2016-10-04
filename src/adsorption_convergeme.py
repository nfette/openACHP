# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 08:21:01 2016

@author: nfette

Adsorption cycle simulation, fast optimization, and parametric study.
"""

from numpy import linspace, zeros, meshgrid, nan, logical_and, unravel_index, array
from scipy.interpolate import PchipInterpolator
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, plot, show, xlabel, ylabel, contour, clabel, title, subplot, xlim, draw
from mpl_toolkits.mplot3d import Axes3D
from adsorption import AdsorptionChiller, AdsorptionChillerControl, AdsorptionChillerSpec, WaterTQ2H
from hw2_1 import KelvinToCelsius as K2C, CelsiusToKelvin as C2K

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

    def __init__(self,ch,q_min=None,q_max=None,refine=False):
        t_cond = ch.ctrl.t_cond
        t_evap = ch.ctrl.t_evap
        t_cool = ch.ctrl.t_cool
        t_exhaust = ch.ctrl.t_exhaust
        
        # The extreme ranges depend on the heat and cold sources.
        margin_up = 1e-3
        margin_down = 1e-5
        if q_min == None:
            self.q_min = ch.f.Q(t_cond, t_exhaust) + margin_down
        else:
            self.q_min = q_min
            
        if q_max == None:
            self.q_max = ch.f.Q(t_evap, t_cool) - margin_up
        else:
            self.q_max = q_max
            
        
        q_range = linspace(self.q_min,self.q_max);
        # We should start desorption at q_max.
        T_d0 = ch.f.T(t_cond, self.q_max);
        T_d1 = ch.f.T(t_cond, self.q_min);
        
        T_a0 = ch.f.T(t_evap, self.q_min);
        T_a1 = ch.f.T(t_evap, self.q_max);
        
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
        # Time functions
        self.ppd = PchipInterpolator(self.q_d[::-1], self.t_d[::-1])
        self.ppa = PchipInterpolator(self.q_a, self.t_a)
        # Compression / decompression is time delta only
        self.pp_comp = PchipInterpolator(q_range,deltat_cmp)
        self.pp_decomp = PchipInterpolator(q_range,deltat_dec)
        # Temperature functions
        self.pp_Td = PchipInterpolator(self.q_d[::-1],T_d[::-1])
        self.pp_Ta = PchipInterpolator(self.q_a,T_a)
        # q as function of T, required to integrate for heat input
        self.pp_q_of_T = PchipInterpolator(T_d,self.q_d)
        self.pp_qintegral = self.pp_q_of_T.antiderivative()
        
        # Stuff to display for user
        td_splined = self.ppd(q_range)
        ta_splined = self.ppa(q_range)
        
        figure(1)
        if not refine:
            plt.cla()
        plot(self.t_d,self.q_d,'ro',label='Desorb')
        plot(self.t_a,self.q_a,'bs',label='Adsorb')
        plot(deltat_cmp,q_range,'pink',label='Compress')
        plot(deltat_dec,q_range,'g',label='Decompress')
        plot(td_splined,q_range,'r')
        plot(ta_splined,q_range,'b')
        xlabel('Integration time, $t$ (s)')
        ylabel('Adsorbed mass ratio, $q$ (kg/kg)')
        title('$T_{{\\rm exhaust}}$ = {:g} K'.format(t_exhaust))
        if not refine:
            plt.legend()
        
    def parametric(self,ch,refine=False):
        """Now do a parametric study. Choose q_low and q_high, then determine
        the time needed in each step, and the other outputs of the system."""
        
        if True:
            Ni = 200
            Nj = 250
            self.q_lows = linspace(self.q_min + 0.0001, self.q_max - 0.001, Ni)
            self.q_highs = linspace(self.q_min + 0.0001, self.q_max - 0.001, Nj)
            #self.q_lows = linspace(0.02, 0.05, Ni)
            #self.q_highs = linspace(0.07, 0.11, Nj)
        else:
            self.q_lows = self.q_d[::2]
            self.q_highs = self.q_a[::2]
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
        # X is low, Y is high
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
        # End of desorption: P_cond, q_low
        T_d1 = self.pp_Td(X)
        # Skip decompression ...
        # End of adsorption: P_evap, q_high
        T_a1 = self.pp_Ta(Y)
        # End of compression: P_cond, q_high
        T_c1 = self.pp_Td(Y)
        
        m_ref = ch.spec.N_beds * (Y - X) * ch.spec.m2 / self.t_cycle
        DeltaH = WaterTQ2H(ch.ctrl.t_evap,1) - WaterTQ2H(ch.ctrl.t_cond,0)
        self.q_ref = m_ref * DeltaH
        # Changed: the Y and T_ terms require an integral.
        dead_mass_capacity = (ch.spec.m2 * ch.spec.c2 + ch.spec.m1 * ch.spec.c1)
        
        Q_dot_in_compress = (ch.spec.N_beds / self.t_cycle) \
            * dead_mass_capacity * (T_c1 - T_a1)
        Q_dot_in_desorb = (m_ref * ch.spec.hads) \
            + (ch.spec.N_beds / self.t_cycle) \
            * (dead_mass_capacity * (T_d1 - T_c1) \
               + ch.spec.m2 * ch.spec.cw * (self.pp_qintegral(T_d1) - self.pp_qintegral(T_c1)))
        Q_dot_in = Q_dot_in_compress + Q_dot_in_desorb
        
        self.q_in = Q_dot_in
        self.cop = self.q_ref / self.q_in
        
        # Do not use this equation, it is wrong:
        #t_cycle = t_desorb + t_adsorb;
        # Desorb step fraction of the total time
        t_fractiond = self.t_desorb / self.t_cycle;
        self.mask=logical_and(self.t_desorb >= 0, self.t_adsorb >= 0)
        mask1[self.mask] = 1;
        masknan[self.mask] = 1;
        # This gives a quick guess at the maximum
        Iflat = (self.q_ref * mask1).argmax()
        I1,I2 = unravel_index(Iflat,self.q_ref.shape)
        
        t_d_opt = self.t_desorb[I1,I2]
        t_e_opt = self.t_decomp[I1,I2]
        t_a_opt = self.t_adsorb[I1,I2]
        t_c_opt = self.t_compress[I1,I2]
        Q_max = self.q_ref[I1,I2]
        COP_max = self.cop[I1,I2]
        t_total = self.t_cycle[I1,I2]
        q_low = self.q_lows[I1-1]
        q_high = self.q_highs[I2+1]
        
        if refine:
            from scipy.interpolate import RectBivariateSpline
            from scipy.optimize import minimize
            spline = RectBivariateSpline(self.q_lows,self.q_highs,self.q_ref)
            def fun(x):
                return -spline(x[0],x[1])
            opt = minimize(fun, array([[q_low,q_high]]))
            q_low,q_high = opt.x
            Q_max = -opt.fun
            t_d_opt = self.ppd(q_low) - self.ppd(q_high)
            t_e_opt = self.pp_decomp(q_low)
            t_a_opt = self.ppa(q_high) - self.ppa(q_low)
            t_c_opt = self.pp_comp(q_high)
            t_total = t_d_opt + t_e_opt + t_a_opt + t_c_opt
            T_d1_opt = self.pp_Td(q_high)
            T_a1_opt = self.pp_Ta(q_low)
            m_ref_opt = 2 * (q_high - q_low) * ch.spec.m2 / t_total
            q_in_opt = (m_ref_opt * ch.spec.hads) \
                + 2 * (ch.spec.m2 * (ch.spec.c2 + q_high * ch.spec.cw) + ch.spec.m1 * ch.spec.c1) * (T_d1_opt - T_a1_opt) / (t_total)
            COP_max = Q_max / q_in_opt
        
        if False:
            figure(2)
            if not refine:
                plt.cla()
            #CS = plt.contourf(X,Y,self.t_cycle)
            #clabel(CS, inline=1, fontsize=10)
            plt.pcolormesh(X,Y,self.t_cycle*masknan,vmin=0,vmax=self.t_cycle.max())
            xlabel('Low mass ratio, $q_{low}$ (kg/kg)')
            ylabel('High mass ratio, $q_{high}$ (kg/kg)')
            title('Cycle time, $t_{cycle}$ (s)')
            
            fig=figure(3)
            ax = fig.gca()
            if not refine:
                ax.cla()
            #CS = plt.contourf(X,Y,t_fractiond*masknan)
            #clabel(CS, inline=1, fontsize=10)
            plt.pcolormesh(X,Y,masknan*t_fractiond,vmin=0,vmax=1)
            xlabel('Low mass ratio, $q_{low}$ (kg/kg)')
            ylabel('High mass ratio, $q_{high}$ (kg/kg)')
            title('Desorption step fraction, $t_{des}/t_{cycle}$')
            
            fig = figure(4)
            if not refine:
                plt.clf()
            """ax = fig.add_subplot(111, projection='3d')
            ax.grid(True)
            ax.plot_wireframe(self.t_desorb*masknan,self.t_adsorb*masknan,X,color='b')
            ax.plot_wireframe(self.t_desorb*masknan,self.t_adsorb*masknan,Y,color='r')
            #plot3([t_d_opt t_d_opt], [t_a_opt t_a_opt], [X(I1,I2) Y(I1,I2)], 'ko-')
            #set(gca,'xlim',[0, max(t_desorb(:))])
            #set(gca,'ylim',[0, max(t_adsorb(:))])
            xlabel('Desorb step time, $\Delta t$ (s)')
            ylabel('Adsorb step time, $\Delta t$ (s)')
            ax.set_zlabel('Mass ratio, $q$ (kg/kg)')"""
            ax = fig.add_subplot(111)
            ax.pcolormesh(X,Y,mask1*self.q_ref)
            xlabel('Low mass ratio, $q_{low}$ (kg/kg)')
            ylabel('High mass ratio, $q_{high}$ (kg/kg)')
            title('Cooling capacity, $Q_{cool}$ (kW)')
            
            fig = figure(5)
            if not refine:
                plt.clf()
            if False:
                ax = fig.add_subplot(111, projection='3d')
                ax.grid(True)
                ax.plot_wireframe(self.t_desorb*masknan, self.t_adsorb*masknan, self.q_ref*masknan, color='b')
                ax.plot([t_d_opt, t_d_opt], [t_a_opt, t_a_opt], [0, Q_max], 'ko-')
                ax.set_zlabel('Cooling capacity, $Q_{cool}$ (kW)')
            else:
                #ax.plot_surface(self.t_desorb*masknan, self.t_adsorb*masknan, self.q_ref*masknan)
                ax = fig.add_subplot(111)
                #ax.contour(self.t_desorb*masknan, self.t_adsorb*masknan, self.q_ref*masknan)
                #ax.tripcolor((self.t_desorb*masknan).flat, (self.t_adsorb*masknan).flat, (self.q_ref*masknan).flat)
                # Uncomment this if you want
                tcf = ax.tricontourf(self.t_desorb[self.mask], self.t_adsorb[self.mask], self.q_ref[self.mask])
                if not refine:
                    fig.colorbar(tcf)
                    ax.set_xlim([0,self.t_desorb.max()])
                    ax.set_ylim([0,self.t_adsorb.max()])
                ax.set_title('Cooling capacity, $Q_{cool}$ (kW)')
            ax.set_xlabel('Desorb step time, $\Delta t$ (s)')
            ax.set_ylabel('Adsorb step time, $\Delta t$ (s)')
        
        t_opts = (t_d_opt, t_e_opt, t_a_opt, t_c_opt)
        return t_opts, Q_max, COP_max, t_total, q_low, q_high

class parstudy1():
    def __init__(self):
        #global t_exhaust T_exhaust_range
        spec = AdsorptionChillerSpec()
        ctrl = AdsorptionChillerControl()
        self.chiller = AdsorptionChiller(spec, ctrl)
        
        Ni = 20;
        #T_exhaust_range = linspace(313,393,Ni);
        self.T_exhaust_range = linspace(313,423,Ni)
        self.tt_d = zeros(Ni) * nan
        self.tt_e = zeros(Ni) * nan
        self.tt_a = zeros(Ni) * nan
        self.tt_c = zeros(Ni) * nan
        self.QQ = zeros(Ni) * nan
        self.CC = zeros(Ni) * nan
        self.tt_total = zeros(Ni)
        
        fig=figure(6)
        plt.clf()
        x = K2C(self.T_exhaust_range)
        self.ax1=subplot(3,1,1)
        self.h = plot(x,self.QQ,'ko')[0]
        self.ax1.relim()
        self.ax1.autoscale_view()
        ylabel('$Q_{\\rm max}$ (kW)')
        self.ax2=subplot(3,1,2,sharex=self.ax1)
        self.h2 = plot(x,self.CC,'ko')[0]
        self.ax2.relim()
        self.ax2.autoscale_view()
        ylabel('COP')
        self.ax3=subplot(3,1,3,sharex=self.ax1)
        self.ax3.stackplot(x,self.tt_d,self.tt_e,self.tt_a,self.tt_c,colors=['r','pink','b','g'])
        self.ax3.set_xlim([x.min(),x.max()])
        self.ax3.set_ylim([0,60])
        self.ax3.set_ylabel('Cycle time (s)')
        self.ax3.set_xlabel('Exhaust temperature, $T$ ($^\circ$C)')
        self.ax3.hold(False)
        self.it = iter(enumerate(self.T_exhaust_range))
    
    def __iter__(self):
        return self
        
    def __next__(self):
        i,t_exhaust = self.it.__next__()
        
        self.chiller.ctrl.t_exhaust = t_exhaust
        print('{} K, '.format(t_exhaust))
        self.something = convergeMe(self.chiller)
        t_opts, Q_max, COP_max, t_total, _, _\
            = self.something.parametric(self.chiller)
        self.tt_d[i],self.tt_e[i],self.tt_a[i],self.tt_c[i] = t_opts
        self.tt_total[i] = t_total
        self.QQ[i] = Q_max
        self.CC[i] = COP_max
        print('{:g} s, {:g} s, {:g} s, {:g} s, {:g} kW\n'.format(*t_opts, Q_max))
    
    def myplot(self):
        self.h.set_ydata(self.QQ)
        self.ax1.relim()
        self.ax1.autoscale_view()
        self.h2.set_ydata(self.CC)
        self.ax2.relim()
        self.ax2.autoscale_view()
        #self.ax3.cla()
        x = K2C(self.T_exhaust_range)
        self.ax3.stackplot(x,self.tt_d,self.tt_e,self.tt_a,self.tt_c,colors=['r','pink','b','g'])
        self.ax3.set_ylim([0,max(self.tt_total)])
        
    def saveplots(self):
        for f in range(1,7):
            fig = figure(f)
            fig.canvas.draw()
            fig.canvas.flush_events()
            plt.savefig('../img/convergeme.{:02d}.{:02d}.png'.format(f,i))

class parstudy2():
    def __init__(self):
        #global t_exhaust T_exhaust_range
        spec = AdsorptionChillerSpec()
        ctrl = AdsorptionChillerControl()
        self.chiller = AdsorptionChiller(spec, ctrl)
        
        Ni = 20;
        #T_exhaust_range = linspace(313,393,Ni);
        self.u_des_range = linspace(0.1,1,Ni)
        self.tt_d = zeros(Ni) * nan
        self.tt_e = zeros(Ni) * nan
        self.tt_a = zeros(Ni) * nan
        self.tt_c = zeros(Ni) * nan
        self.QQ = zeros(Ni) * nan
        self.CC = zeros(Ni) * nan
        self.tt_total = zeros(Ni)
        self.qlow = zeros(Ni) * nan
        self.qhigh = zeros(Ni) * nan
        
        fig=figure(6)
        plt.clf()
        X = self.u_des_range
        self.ax1=subplot(3,1,1)
        self.h = plot(X,self.QQ,'ko')[0]
        self.ax1.relim()
        self.ax1.autoscale_view()
        ylabel('$Q_{\\rm max}$ (kW)')
        self.ax2=subplot(3,1,2,sharex=self.ax1)
        self.h2 = plot(X,self.CC,'ko')[0]
        self.ax2.relim()
        self.ax2.autoscale_view()
        ylabel('COP')
        self.ax3=subplot(3,1,3,sharex=self.ax1)
        self.ax3.stackplot(X,self.tt_d,self.tt_e,self.tt_a,self.tt_c,colors=['r','pink','b','g'])
        self.ax3.set_xlim([X.min(),X.max()])
        self.ax3.set_ylim([0,60])
        self.ax3.set_ylabel('Cycle time (s)')
        self.ax3.set_xlabel(r'Desorber overall HX coefficient, $U$ (kW$\cdot$ m$^{-2}\cdot$ K$^{-1}$)')
        self.ax3.hold(False)
        self.it = iter(enumerate(X))
    
    def __iter__(self):
        return self
        
    def __next__(self):
        i,x = self.it.__next__()
        
        self.chiller.spec.u_des = x
        print('{} kW/m^2-K, '.format(x))
        # Course result and plots
        self.something = convergeMe(self.chiller)
        t_opts, Q_max, COP_max, t_total, q_low, q_high\
            = self.something.parametric(self.chiller,False)
        print('{:g} s, {:g} s, {:g} s, {:g} s, {:g} kW, {:g}, {:g}\n'.format(*t_opts, Q_max, q_low, q_high))
        self.tt_d[i],self.tt_e[i],self.tt_a[i],self.tt_c[i] = t_opts
        self.tt_total[i] = t_total
        self.QQ[i] = Q_max
        self.CC[i] = COP_max
        self.qlow[i] = q_low
        self.qhigh[i] = q_high
        
        self.h.set_ydata(self.QQ)
        self.ax1.relim()
        self.ax1.autoscale_view()
        self.h2.set_ydata(self.CC)
        self.ax2.relim()
        self.ax2.autoscale_view()
        #self.ax3.cla()
        X = self.u_des_range
        self.ax3.stackplot(X,self.tt_d,self.tt_e,self.tt_a,self.tt_c,colors=['r','pink','b','g'])
        self.ax3.set_ylim([0,max(self.tt_total)])
        for f in range(1,7):
            fig = figure(f)
            fig.canvas.draw()
            fig.canvas.flush_events()
            plt.savefig('../img/convergeme2.{:02d}.{:02d}.png'.format(f,i))
            
if __name__ == "__main__":
    #input('Hit enter to proceed')
    #ps = parstudy1()
    import pickle
    try:
        with open('mypickle.pkl','rb') as f:
            ps = pickle.load(f)
    except:
        ps = parstudy2()
        for p in ps:
            pass
        with open('mypickle.pkl','wb') as f:
            pickle.dump(ps,f)
    
    ch=ps.chiller
    total_adsorbent_mass = ch.spec.m2 * ch.spec.N_beds
    SCP = ps.QQ / total_adsorbent_mass
    DeltaH = WaterTQ2H(ch.ctrl.t_evap,1) - WaterTQ2H(ch.ctrl.t_cond,0)
    Deltaq = ps.qhigh - ps.qlow
    SCPapprox = DeltaH * Deltaq / ps.tt_total
#    fig,ax=plt.subplots()
#    mylabel=r'Desorber overall HX coefficient, $U$ (kW$\cdot$ m$^{-2}\cdot$ K$^{-1}$)'
#    ax.set_xlabel(mylabel)
#    ax.set_ylabel('SCP [kW/kg]')
#    ax.grid()
#    ax.plot(ps.u_des_range, SCP)
    ax=plt.twinx(ps.ax1)
    ax.plot(ps.u_des_range,SCP)
    ax.axis('tight')
    ps.ax1.axis('tight')
    ax.set_ylabel('SCP [kW/kg]')
    
    ps.ax1.grid()
    ps.ax2.grid()
    ps.ax3.grid()
    
    ps.h.set_marker(None); ps.h.set_linestyle('-')
    ps.h2.set_marker(None); ps.h2.set_linestyle('-')
    plt.gcf().tight_layout()
    plt.draw()