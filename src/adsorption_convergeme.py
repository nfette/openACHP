# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 08:21:01 2016

@author: nfette

Adsorption cycle simulation, fast optimization, and parametric study.
"""

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
def convergeMe(t_cond, t_evap, t_cool, t_exhaust):
    
    # The extreme ranges depend on the heat and cold sources.
    q_min = freundlich(t_cond, t_exhaust);
    q_max = freundlich(t_evap, t_cool);
    # We should start desorption at q_max.
    T_d0 = Tfreundlich(t_cond, q_max-0.001);
    T_d1 = Tfreundlich(t_cond, q_min+0.0001);
    
    T_a0 = Tfreundlich(t_evap, q_min+0.0001);
    T_a1 = Tfreundlich(t_evap, q_max-0.001);
    
    #td = linspace(0,200,100);
    #[ty, qy] = desorption9(td, T_d0);
    T_d = linspace(T_d0,T_d1);
    [t_d, q_d] = desorptionFlip(T_d);
    
    #ta = linspace(0,1000,100);
    #[tad, qad] = adsorption9(ta, T_a0);
    T_a = linspace(T_a0, T_a1);
    [t_a, q_a] = adsorptionFlip(T_a);
    
    # Spline it.
    ppd = interp1(q_d,t_d,'pchip','pp');
    ppa = interp1(q_a,t_a,'pchip','pp');
    q_range = linspace(q_min+0.0001,q_max-0.01);
    td_splined = ppval(ppd,q_range);
    ta_splined = ppval(ppa,q_range);
    
    figure(1)
    plot(t_d,q_d,'ro',t_a,q_a,'bo',td_splined,q_range,'r',ta_splined,q_range,'b')
    xlabel('Integration time, $$t$$ (s)')
    ylabel('Adsorbed mass ratio, $$q$$ (kg/kg)')
    title(sprintf('$$T_{\\rm exhaust}$$ = %g K',t_exhaust))
    
    # Now do a parametric study. Choose q_low and q_high, then determine
    # the time needed in each step, and the other outputs of the system.
    
    Ni = 100;
    Nj = 100;
    q_lows = linspace(q_min + 0.0001, q_max - 0.001, Ni);
    q_highs = linspace(q_min + 0.0001, q_max - 0.001, Nj);
    q_lows = q_d(1:5:end);
    q_highs = q_a(1:5:end);
    Ni = length(q_lows);
    Nj = length(q_highs);
    t_desorb = nan(Ni, Nj);
    t_decomp = nan(Ni, Nj);
    t_adsorb = nan(Ni, Nj);
    t_compress = nan(Ni, Nj);
    
    t_cycle = nan(Ni, Nj);
    q_ref = nan(Ni, Nj);
    q_in = nan(Ni, Nj);
    cop = nan(Ni, Nj);
    mask = nan(Ni, Nj);
    for i = 1:Ni
    q_low = q_lows(i);
    for j = 1:Nj
        q_high = q_highs(j);
        if q_high < q_low
            #continue
        end
        t_desorb(i,j) = diff(ppval(ppd,[q_high,q_low]));
        t_decomp(i,j) = decompress(q_low);
        t_adsorb(i,j) = diff(ppval(ppa,[q_low,q_high]));
        t_compress(i,j) = compress(q_high);
        t_cycle(i,j) = t_desorb(i,j) + t_compress(i,j) + t_adsorb(i,j) + t_decomp(i,j);
        [q_ref(i,j), q_in(i,j), cop(i,j)] = afterSolve(q_low,q_high,t_cycle(i,j));
    end
    end
    #t_cycle = t_desorb + t_adsorb;
    t_fractiond = t_desorb ./ t_cycle;
    mask(t_desorb >= 0 & t_adsorb >= 0) = 1;
    [M,I] = max(q_ref .* mask);
    [M2,I2] = max(M);
    I1 = I(I2);
    t_d_opt = t_desorb(I1,I2);
    t_e_opt = t_decomp(I1,I2);
    t_a_opt = t_adsorb(I1,I2);
    t_c_opt = t_compress(I1,I2);
    Q_max = M2;
    COP_max = cop(I1,I2);
    
    figure(2)
    [X,Y] = ndgrid(q_lows,q_highs);
    #surf(X,Y,t_desorb,'LineStyle','none');hold on
    #surf(X,Y,t_adsorb,'LineStyle','none')
    contour(X,Y,t_cycle,'ShowText','on')
    xlabel('Low mass ratio, $$q_low$$ (kg/kg)')
    ylabel('High mass ratio, $$q_high$$ (kg/kg)')
    zlabel('Cycle time, $$t_cycle$$ (s)')
    #h = rotate3d;
    #set(h,'ActionPostCallback',@align_axislabels)
    #set(gca,'DataAspectRatio', get(gca,'DataAspectRatio'))
    
    figure(3)
    contour(X,Y,t_fractiond,'ShowText','on')
    xlabel('Low mass ratio, $$q_low$$ (kg/kg)')
    ylabel('High mass ratio, $$q_high$$ (kg/kg)')
    zlabel('Cycle time, $$t_cycle$$ (s)')
    
    figure(4)
    plot3(t_desorb, t_adsorb, X, 'b-'); hold on
    plot3(t_desorb, t_adsorb, Y, 'r-')
    plot3(t_desorb', t_adsorb', X', 'b-')
    plot3(t_desorb', t_adsorb', Y', 'r-')
    plot3([t_d_opt t_d_opt], [t_a_opt t_a_opt], [X(I1,I2) Y(I1,I2)], 'ko-')
    hold off
    grid on
    set(gca,'xlim',[0, max(t_desorb(:))])
    set(gca,'ylim',[0, max(t_adsorb(:))])
    xlabel('Desorb step time, $$\Delta t$$ (s)')
    ylabel('Adsorb step time, $$\Delta t$$ (s)')
    zlabel('Mass ratio, $$q$$ (kg/kg)')
    h = rotate3d;
    set(h,'ActionPostCallback',@align_axislabels)
    set(gca,'DataAspectRatio', get(gca,'DataAspectRatio'))
    align_axislabels([],h)
    
    figure(5)
    plot3(t_desorb.*mask, t_adsorb.*mask, q_ref.*mask, 'b-'); hold on
    plot3((t_desorb.*mask)', (t_adsorb.*mask)', (q_ref.*mask)', 'r-')
    plot3([t_d_opt t_d_opt], [t_a_opt t_a_opt], [0 Q_max], 'ko-')
    hold off
    #plot3(t_desorb, t_adsorb, q_in, 'b-')
    #plot3(t_desorb', t_adsorb', q_in', 'r-')
    #set(gca,'xlim',[0, max(t_desorb(:))])
    #set(gca,'ylim',[0, max(t_adsorb(:))])
    grid on
    #set(gca,'zlim',[0, 12])
    xlabel('Desorb step time, $$\Delta t$$ (s)')
    ylabel('Adsorb step time, $$\Delta t$$ (s)')
    zlabel('Cooling capacity, $$Q_{cool}$$ (kW)')
    h = rotate3d;
    set(h,'ActionPostCallback',@align_axislabels)
    set(gca,'DataAspectRatio', get(gca,'DataAspectRatio'))
    align_axislabels([],h)
    
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
    global t_exhaust T_exhaust_range
    
    Ni = 21;
    T_exhaust_range = linspace(313,393,Ni);
    tt_d = nan(Ni,1);
    tt_e = nan(Ni,1);
    tt_a = nan(Ni,1);
    tt_c = nan(Ni,1);
    QQ = nan(Ni,1);
    CC = nan(Ni,1);
    
    figure(6)
    subplot(2,1,1)
    h = plot(T_exhaust_range-273,QQ,'ko');
    xlim([40,120])
    ylabel('$$Q_{\rm max}$$ (kW)')
    subplot(2,1,2)
    h2 = plot(T_exhaust_range-273,CC,'ko');
    xlim([40,120])
    ylabel('COP')
    xlabel('Exhaust temperature, $$T$$ ($$^\circ$$C)')
    
    
    ##
    for i = 1:Ni
        t_exhaust = T_exhaust_range(i);
        fprintf('%g K, ',t_exhaust)
        [t_d_opt, t_e_opt, t_a_opt, t_c_opt, Q_max, COP_max] = convergeMe();
        tt_d(i) = t_d_opt;
        tt_e(i) = t_e_opt;
        tt_a(i) = t_a_opt;
        tt_c(i) = t_c_opt;
        QQ(i) = Q_max;
        CC(i) = COP_max;
        fprintf('%g s, %g s, %g s, %g s, %g kW\n',t_d_opt, t_e_opt, t_a_opt, t_c_opt, Q_max)
        set(h, 'YData', QQ)
        set(h2, 'YData', CC)
        for f = 1:5
            figure(f)
            drawnow
        end
    end
    
    %%
