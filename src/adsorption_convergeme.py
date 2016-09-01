# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 08:21:01 2016

@author: nfette

Adsorption cycle simulation, fast optimization, and parametric study.
"""

from numpy import linspace, zeros, diff, meshgrid, nan, logical_and, unravel_index
from scipy.interpolate import PchipInterpolator
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, plot, show, xlabel, ylabel, contour, clabel, title, subplot, xlim, draw
from mpl_toolkits.mplot3d import Axes3D
from adsorption import * 

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
def convergeMe(self):
    t_cond = self.ctrl.t_cond
    t_evap = self.ctrl.t_evap
    t_cool = self.ctrl.t_cool
    t_exhaust = self.ctrl.t_exhaust
    
    # The extreme ranges depend on the heat and cold sources.
    q_min = self.f.Q(t_cond, t_exhaust);
    q_max = self.f.Q(t_evap, t_cool);
    # We should start desorption at q_max.
    T_d0 = self.f.T(t_cond, q_max-0.001);
    T_d1 = self.f.T(t_cond, q_min+0.0001);
    
    T_a0 = self.f.T(t_evap, q_min+0.0001);
    T_a1 = self.f.T(t_evap, q_max-0.001);
    
    #td = linspace(0,200,100);
    #[ty, qy] = desorption9(td, T_d0);
    T_d = linspace(T_d0,T_d1);
    [t_d, q_d] = self.desorptionFlip(T_d);
    
    #ta = linspace(0,1000,100);
    #[tad, qad] = adsorption9(ta, T_a0);
    T_a = linspace(T_a0, T_a1);
    [t_a, q_a] = self.adsorptionFlip(T_a);
    
    # Spline it.
    ppd = PchipInterpolator(q_d[::-1], t_d[::-1])
    ppa = PchipInterpolator(q_a, t_a)
    q_range = linspace(q_min+0.0001,q_max-0.01);
    td_splined = ppd(q_range)
    ta_splined = ppa(q_range);
    
    figure(1)
    plot(t_d,q_d,'ro',t_a,q_a,'bo',td_splined,q_range,'r',ta_splined,q_range,'b')
    xlabel('Integration time, $t$ (s)')
    ylabel('Adsorbed mass ratio, $q$ (kg/kg)')
    title('$T_{{\\rm exhaust}}$ = {:g} K'.format(t_exhaust))
    
    # Now do a parametric study. Choose q_low and q_high, then determine
    # the time needed in each step, and the other outputs of the system.
    
    if False:
        Ni = 100
        Nj = 100
        q_lows = linspace(q_min + 0.0001, q_max - 0.001, Ni)
        q_highs = linspace(q_min + 0.0001, q_max - 0.001, Nj)
    else:
        q_lows = q_d[::5]
        q_highs = q_a[::5]
        Ni = len(q_lows)
        Nj = len(q_highs)
    t_desorb = zeros((Ni, Nj))
    t_decomp = zeros((Ni, Nj))
    t_adsorb = zeros((Ni, Nj))
    t_compress = zeros((Ni, Nj))
    
    t_cycle = zeros((Ni, Nj))
    q_ref = zeros((Ni, Nj))
    q_in = zeros((Ni, Nj))
    cop = zeros((Ni, Nj))
    mask = zeros((Ni, Nj))
    mask.fill(nan)
    for i, q_low in enumerate(q_lows):
        for j, q_high in enumerate(q_highs):
            if q_high < q_low:
                #continue
                pass
            t_desorb[i,j] = ppd(q_low) - ppd(q_high) 
            t_decomp[i,j] = self.decompress(q_low)
            t_adsorb[i,j] = ppa(q_high) - ppa(q_low)
            t_compress[i,j] = self.compress(q_high)
            t_cycle[i,j] = t_desorb[i,j] + t_compress[i,j] + t_adsorb[i,j] + t_decomp[i,j]
            q_ref[i,j], q_in[i,j], cop[i,j], _ = self.afterSolve(q_low,q_high,t_cycle[i,j])
    
    # Do not use this equation, it is wrong:
    #t_cycle = t_desorb + t_adsorb;
    # Desorb step fraction of the total time
    t_fractiond = t_desorb / t_cycle;
    mask[logical_and(t_desorb >= 0, t_adsorb >= 0)] = 1;
    #[M,I] = max(q_ref .* mask);
    #[M2,I2] = max(M);
    #I1 = I(I2);
    Iflat = (q_ref * mask).argmax()
    I1,I2 = unravel_index(Iflat,q_ref.shape)    
    t_d_opt = t_desorb[I1,I2];
    t_e_opt = t_decomp[I1,I2];
    t_a_opt = t_adsorb[I1,I2];
    t_c_opt = t_compress[I1,I2];
    Q_max = q_ref[I1,I2];
    COP_max = cop[I1,I2];
    
    figure(2)
    [X,Y] = meshgrid(q_lows,q_highs,indexing='ij')
    CS = contour(X,Y,t_cycle)
    clabel(CS, inline=1, fontsize=10)    
    xlabel('Low mass ratio, $q_low$ (kg/kg)')
    ylabel('High mass ratio, $q_high$ (kg/kg)')
    title('Cycle time, $t_cycle$ (s)')
    
    figure(3)
    CS = contour(X,Y,t_fractiond)
    clabel(CS, inline=1, fontsize=10)
    xlabel('Low mass ratio, $q_low$ (kg/kg)')
    ylabel('High mass ratio, $q_high$ (kg/kg)')
    title('Desorption step fraction, $t_des/t_cycle$')
    
    fig = figure(4)
    ax = fig.add_subplot(111, projection='3d')
    ax.grid(True)
    #plot3(t_desorb, t_adsorb, X, 'b-'); hold on
    ax.plot_wireframe(t_desorb,t_adsorb,X,color='b')
    ax.plot_wireframe(t_desorb,t_adsorb,Y,color='r')
    #plot3(t_desorb, t_adsorb, Y, 'r-')
    #plot3(t_desorb.T, t_adsorb.T, X.T, 'b-')
    #plot3(t_desorb.T, t_adsorb.T, Y.T, 'r-')
    #plot3([t_d_opt t_d_opt], [t_a_opt t_a_opt], [X(I1,I2) Y(I1,I2)], 'ko-')
    #set(gca,'xlim',[0, max(t_desorb(:))])
    #set(gca,'ylim',[0, max(t_adsorb(:))])
    xlabel('Desorb step time, $\Delta t$ (s)')
    ylabel('Adsorb step time, $\Delta t$ (s)')
    ax.set_zlabel('Mass ratio, $q$ (kg/kg)')
    
    fig = figure(5)
    ax = fig.add_subplot(111, projection='3d')
    ax.grid(True)
    ax.plot_wireframe(t_desorb*mask, t_adsorb*mask, q_ref*mask, color='b')
    #plot3((t_desorb.*mask).T, (t_adsorb.*mask).T, (q_ref.*mask).T, 'r-')
    ax.plot([t_d_opt, t_d_opt], [t_a_opt, t_a_opt], [0, Q_max], 'ko-')
    #plot3(t_desorb, t_adsorb, q_in, 'b-')
    #plot3(t_desorb', t_adsorb', q_in', 'r-')
    #set(gca,'xlim',[0, max(t_desorb(:))])
    #set(gca,'ylim',[0, max(t_adsorb(:))])    
    #set(gca,'zlim',[0, 12])
    xlabel('Desorb step time, $\Delta t$ (s)')
    ylabel('Adsorb step time, $\Delta t$ (s)')
    ax.set_zlabel('Cooling capacity, $Q_{cool}$ (kW)')
    #h = rotate3d;
    #set(h,'ActionPostCallback',@align_axislabels)
    #set(gca,'DataAspectRatio', get(gca,'DataAspectRatio'))
    #align_axislabels([],h)
    
    # figure()
    # plot3(t_desorb, t_adsorb, cop, 'b-'); hold on
    # plot3(t_desorb', t_adsorb', cop', 'r-')
    # #set(gca,'xlim',[0, max(t_desorb(:))])
    # #set(gca,'ylim',[0, max(t_adsorb(:))])
    # set(gca,'zlim',[-1, 1])
    # xlabel('Desorb step time')
    # ylabel('Adsorb step time')
    # zlabel('COP')

    return [t_d_opt, t_e_opt, t_a_opt, t_c_opt, Q_max, COP_max]

if __name__ == "__main__":
    #global t_exhaust T_exhaust_range
    spec = AdsorptionChillerSpec()
    ctrl = AdsorptionChillerControl()
    chiller = AdsorptionChiller(spec, ctrl)
    
    Ni = 5;
    T_exhaust_range = linspace(313,393,Ni);
    tt_d = zeros(Ni)
    tt_e = zeros(Ni)
    tt_a = zeros(Ni)
    tt_c = zeros(Ni)
    QQ = zeros(Ni)
    CC = zeros(Ni)
    
    figure(6)
    subplot(2,1,1)
    h = plot(T_exhaust_range-273,QQ,'ko')[0]
    xlim([40,120])
    ylabel('$Q_{\\rm max}$ (kW)')
    subplot(2,1,2)
    h2 = plot(T_exhaust_range-273,CC,'ko')[0]
    xlim([40,120])
    ylabel('COP')
    xlabel('Exhaust temperature, $T$ ($^\circ$C)')
    
    
    ##
    for i,t_exhaust in enumerate(T_exhaust_range):
        ctrl.t_exhaust = t_exhaust
        print('{} K, '.format(t_exhaust))
        t_d_opt, t_e_opt, t_a_opt, t_c_opt, Q_max, COP_max = convergeMe(chiller)
        tt_d[i] = t_d_opt
        tt_e[i] = t_e_opt
        tt_a[i] = t_a_opt
        tt_c[i] = t_c_opt
        QQ[i] = Q_max
        CC[i] = COP_max
        print('{:g} s, {:g} s, {:g} s, {:g} s, {:g} kW\n'.format(t_d_opt, t_e_opt, t_a_opt, t_c_opt, Q_max))
        h.set_ydata(QQ)
        h2.set_ydata(CC)
        for f in range(1,6):
            figure(f)
            draw()
        show()
