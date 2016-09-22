# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 19:15:24 2016

@author: nfette

Adsorption cycle simulation (script only).

Translated from Gupta's original MATLAB solution. This is a parametric study
wrt heat input stream temperature and half cycle time (adsorption and
desorption steps are given the same dwell time). Key features:
* Integrates the ODE for adsorbed mass (& temperature) over given time interval
* First, determine the temperature at start of desorption, T_d0, when the
  system has reached a steady cycle.
* Then compute the states, refrigerant mass flow, etc, and output Q and COP.

"""
from __future__ import print_function
from numpy import linspace, logspace, nan, zeros, array, meshgrid
from scipy.optimize import fsolve
from matplotlib.pyplot import plot, figure, xlabel, ylabel, draw
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from adsorption import *

def main():
    spec = AdsorptionChillerSpec()
    ctrl = AdsorptionChillerControl()
    chiller = AdsorptionChiller(spec, ctrl)
    
    Ni, Nj = 2, 3
    T_exhaust_range = linspace(313,393,Ni)
    end_t_range = logspace(2,3,Nj)
    #end_t_range = [240];
    
    #Ni = length(T_exhaust_range);
    #Nj = length(end_t_range);
    q_ref = zeros((Ni, Nj))
    q_in = zeros((Ni, Nj))
    cop = zeros((Ni, Nj))
    T = []
    ##
    T_d0 = 311
    #T_d0 = 325.6137
    
    for i, t_exhaust in enumerate(T_exhaust_range):
        ctrl.t_exhaust = t_exhaust
        #for j = Nj:-1:1
        for j, end_t in enumerate(end_t_range):
            ctrl.end_t = end_t
            t = linspace(ctrl.start_t,ctrl.end_t,endpoint=True)
            #[T_d0,fval,exitflag,output,jacobian] = fsolve(@(T)(loopOnce(T,t)),T_d0);
            x,infodict,ier,mesg = fsolve(chiller.loopOnce,T_d0,args=(t,t),full_output=True)
            print(mesg)
            T_d0 = x
    
            if 'qad' in vars():
                del qad
            print('T={}, '.format(t_exhaust))
            fig3=figure(3)
            fig3.clear()
            xlabel('Time(sec)')
            ylabel('Temperature($^\circ$C)')
            ax3=fig3.gca()
            ax3.cla()
            #hold on
            fig4 = figure(4) # cla
            fig4.clear()
            xlabel('$T$ / $^\circ$C')
            ylabel('$q$ / (kg/kg)')
            ax4=fig4.gca()
            ax4.cla()
            ax4.grid(True)
            #hold on
    
            #       for t_evap = [283:1:295]
            #       for m_water = [0.1:0.01:0.6]
            
            dta = 1
            myt = 0
            for k in range(5):
            #while dta>=1:
            #while dta>=1e-3:
                dTa, dt, q_d1, q_a1 = chiller.loopOnce(T_d0, t, t, ax3, ax4, [fig3,fig4], myt)
                myt += dt
                """
                ty, qy = chiller.desorption9(t, T_d0)
                fig=figure(3)
                plot(myt + t,ty-273.15,'r')
                draw()
                fig.canvas.flush_events()
                figure(4);
                if 'qad' in vars():
                    # There was a compression step
                    # Continue plot from previous adsorption step
                    plot(K2C(array([tad[-1],ty[0]])),[qad[-1],qy[0]],'o-');
                plot(ty-273.15,qy,'r'); draw()
                myt = myt + end_t
                T_d1 = ty[-1]
                q_d1 = chiller.f.Q(ctrl.t_cond, T_d1)
                P_a0 = wsr4t2p(ctrl.t_evap) / ((q_d1 / spec.xo) ** spec.n)
                T_a0 = wsr4p2t(P_a0)
                q_a0 = q_d1
    
                tad, qad = chiller.adsorption9(t, T_a0)
                T_a1 = tad[-1]
                q_a1 = chiller.f.Q(ctrl.t_evap, T_a1);
                P_d0 = wsr4t2p(ctrl.t_cond) / ((q_a1 / spec.xo) ** spec.n)
                T_d0_new = wsr4p2t(P_d0)
    
                figure(3); plot(myt + t,K2C(tad),'b'); draw()
                figure(4); plot(K2C(array([T_d1,T_a0])),[q_d1,q_a0],'o-')
                plot(K2C(tad),qad,'b'); draw()
                myt = myt + end_t
    
                dta = abs(T_d0 - T_d0_new)
                T_d0 = T_d0_new
                fig.canvas.flush_events()
                """
    
            x_dil=q_d1
            x_conc=q_a1
            q_ref_this, q_in_this, cop_this, m_ref = chiller.afterSolve(x_dil, x_conc, 0.5*end_t)
            pe = wsr4t2p(ctrl.t_evap)
    
            q_ref[i,j] = q_ref_this
            q_in[i,j] = q_in_this
            cop[i,j] = cop_this
            print('{}, {}\n'.format(cop[i,j],q_ref[i,j]))
            #raw_input("Please press [Enter] to proceed")
            #plt.draw()
            #plt.show()
            
    T = K2C(T_exhaust_range)
    
    # <codecell>
    fig=figure(7)
    #set(gca,'dataaspectratio',[1, 1, 0.5],'projection','perspective','box','on')
    ax = fig.add_subplot(111,projection='3d')
    
    [X,Y]=meshgrid(T,end_t_range,indexing='ij');
    surf=ax.plot_surface(X,Y,q_ref, cmap=cm.viridis, rstride=1, cstride=1)
    xlabel('Hot Water Temperature, $T$ ($^\circ$C)')
    ylabel('Half cycle time, $t/2$')
    ax.set_zlabel('Cooling capacity, kW')
    fig.colorbar(surf, shrink=0.5, aspect=5)
    # Not necessary in matplotlib
    #set(h,'ActionPostCallback',@align_axislabels)
    #set(gca,'DataAspectRatio', get(gca,'DataAspectRatio'))
    
    figure(1)
    plot(T ,q_ref,'k-'),xlabel('Hot Water Temperature ($^\circ$C)'),ylabel('Refrigeration Capacity(kW)')
    figure(2)
    plot(T ,cop,'k-'),xlabel('Hot Water Temperature ($^\circ$C)'),ylabel('COP')
    
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    #plotyy(T,q_ref,T,cop);
    ax1.plot(T,q_ref,'o:')
    ax2.plot(T,cop,'s:')
    ax1.set_xlabel('Hot Water Temperature ($^\circ$C)')
    ax1.set_ylabel('Refrigeration Capacity (kW)') # left y-axis
    ax2.set_ylabel('COP') # right y-axis
    
    ax1.set_ylim([0,21])
    #ax1.YTickMode='auto';
    ax2.set_ylim([0,1.4])
    #hAx(2).YTickMode='auto';
    
    plt.show()

    
if __name__ == "__main__":
    main()
    