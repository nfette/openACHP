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
from numpy import linspace, logspace
from scipy import fsolve
from matplotlib.pyplot import plot


i = 1;

Ni = 2;
Nj = 2;
T_exhaust_range = linspace(313,393,Ni);
end_t_range = logspace(1,2,Nj);
#end_t_range = [240];

#Ni = length(T_exhaust_range);
#Nj = length(end_t_range);
q_ref = nan(Ni, Nj);
q_in = nan(Ni, Nj);
cop = nan(Ni, Nj);
##
T_d0 = 311;
T_d0 = 325.6137;

for i = 1:Ni
    t_exhaust = T_exhaust_range(i);
    for j = Nj:-1:1
        end_t = end_t_range(j);
        t = start_t:1:end_t;
        [T_d0,fval,exitflag,output,jacobian] = fsolve(@(T)(loopOnce(T,t)),T_d0);
        exitflag
        output
        del qad
        fprintf('T=%f, ',t_exhaust)
        figure(3); cla
        xlabel('Time(sec)')
        ylabel('Temperature($$^\circ$$C)');
        hold on
        figure(4); cla
        xlabel('$$T$$ / $$^\circ$$C')
        ylabel('$$q$$ / (kg/kg)');
        grid on
        hold on

        #       for t_evap = [283:1:295]
        #       for m_water = [0.1:0.01:0.6]
        
        dta = 1;
        myt = 0;
        #while dta>=1
        while dta>=1e-3
            [ty, qy] = desorption9(t, T_d0);
            figure(3); plot(myt + t,ty-273.15,'r'); drawnow
            figure(4);
            if exist('qad','var')
                plot([tad(end),ty(1)]-273.15,[qad(end),qy(1)],'o-');
            end
            plot(ty-273.15,qy,'r'); drawnow
            myt = myt + end_t;
            T_d1 = ty(end);
            q_d1 = freundlich(t_cond, T_d1);
            P_a0 = wsr4t2p(t_evap) / ((q_d1 / xo) ^ n);
            T_a0 = wsr4p2t(P_a0);
            q_a0 = q_d1;

            [tad, qad] = adsorption9(t, T_a0);
            T_a1 = tad(end);
            q_a1 = freundlich(t_evap, T_a1);
            P_d0 = wsr4t2p(t_cond) / ((q_a1 / xo) ^ n);
            T_d0_new = wsr4p2t(P_d0);

            figure(3); plot(myt + t,tad-273.15,'b'); drawnow
            figure(4); plot([T_d1,T_a0]-273.15,[q_d1,q_a0],'o-');
            plot(tad-273.15,qad,'b'); drawnow
            myt = myt + end_t;

            dta = abs(T_d0 - T_d0_new);
            T_d0 = T_d0_new;
        end
        plot(T_a1-273.15,q_a1,'o')

        x_dil=q_d1;
        x_conc=q_a1;
        m_ref = ((x_conc - x_dil)*m2)/(end_t);
        pe = wsr4t2p(t_evap);

        q_ref(i,j) = m_ref*(WaterTQ2H(t_evap,1) - WaterTQ2H(t_cond,0));
        q_in(i,j) = (m_ref*hads) + ((m2*(c2+x_conc*cw) +m1*c1) *(T_d1 - T_a1)/end_t);
        if q_ref(i,j)<0
            q_ref(i,j)=0;
        end
        cop(i,j)=q_ref(i,j)/q_in(i,j);
        fprintf('%f, %f\n',cop(i,j),q_ref(i,j))
        T(i) = t_exhaust - 273;
        
        #             TF(i) = (T(i) * (9/5)) + 32;
        #             m(i) = m_water;
        #           Te(i) = t_evap - 273;
        #             Tc(i) = t_cool - 273;
        #             time(i) = end_t*2;
        
    end
end

##
figure(7)
set(gca,'dataaspectratio',[1 1 0.5],'projection','perspective','box','on')

[X,Y]=ndgrid(T_exhaust_range,end_t_range);
mesh(X,Y,q_ref)
xlabel('Hot Water Temperature, $$T$$ ($$^\circ$$C)')
ylabel('Half cycle time, $$t/2$$')
zlabel('Cooling capacity, kW')
h = rotate3d;
set(h,'ActionPostCallback',@align_axislabels)
set(gca,'DataAspectRatio', get(gca,'DataAspectRatio'))

figure(1)
plot(T ,q_ref,'k-'),xlabel('Hot Water Temperature ($$^\circ$$C)'),ylabel('Refrigeration Capacity(kW)')
figure(2)
plot(T ,cop,'k-'),xlabel('Hot Water Temperature ($$^\circ$$C)'),ylabel('COP')

figure(5)
[hAx,hLine1,hLine2] = plotyy(T,q_ref,T,cop);
xlabel('Hot Water Temperature ($$^\circ$$C)')
ylabel(hAx(1),'Refrigeration Capacity (kW)') # left y-axis
ylabel(hAx(2),'COP') # right y-axis
hLine1.LineStyle = ':';
hLine2.LineStyle = ':';
hLine1.Marker = 'o';
hLine2.Marker = 's';
hAx(1).YLim=[0,21];
hAx(1).YTickMode='auto';
hAx(2).YLim=[0,1.4];
hAx(2).YTickMode='auto';
