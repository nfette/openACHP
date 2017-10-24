import numpy
import CoolProp.CoolProp as CP
import CoolProp.HumidAirProp as HAP
import pandas

h_ref_water = CP.PropsSI('H','T',273.15,'Q',0,'water')

def parallel_evaporation(
    temperature_interface = 300,
    temperature_far_in = 290,
    relhum_far_in = 0.2,
    dryair_mass_flow_in = 0.1,
    water_flow_in = 0.07,
    a_cell = 2.,
    pressure = 1e5,
    cells = 20,
    supersaturated_method=None):
    """For supersaturated method, consider 'hack'."""

    nodal_temperature = []
    nodal_relhum = []
    nodal_humrat = []
    nodal_water_flow = []
    nodal_water_temperature = []
    nodal_interface_humrat = []
    relhum_interface = 1

    nodal_temperature.append(temperature_far_in)
    nodal_relhum.append(relhum_far_in)
    humrat_far_in = HAP.HAPropsSI('HumRat','P',pressure,'T',temperature_far_in,'RH',relhum_far_in)
    nodal_humrat.append(humrat_far_in)
    nodal_water_flow.append(water_flow_in)
    nodal_water_temperature.append(temperature_interface)
    nodal_interface_humrat.append(HAP.HAPropsSI(
        'HumRat','P',pressure,'T',temperature_interface,'RH',relhum_interface))

    for i in range(cells):
        temperature_far = temperature_far_in
        relhum_far = relhum_far_in

        # Convert to dimensions we wish to average
        enthalpy_interface = HAP.HAPropsSI('Hha','P',pressure,'T',temperature_interface,'RH',relhum_interface)
        humrat_interface = HAP.HAPropsSI('HumRat','P',pressure,'T',temperature_interface,'RH',relhum_interface)
        mole_fraction_interface = HAP.HAPropsSI('Y','P',pressure,'T',temperature_interface,'RH',relhum_interface)
        mixture_density_interface = 1 / HAP.HAPropsSI('Vha','P',pressure,'T',temperature_interface,'RH',relhum_interface)
        vapor_density_interface = mole_fraction_interface * mixture_density_interface

        enthalpy_far = HAP.HAPropsSI('Hha','P',pressure,'T',temperature_far,'RH',relhum_far)
        humrat_far = HAP.HAPropsSI('HumRat','P',pressure,'T',temperature_far,'RH',relhum_far)
        mole_fraction_far = HAP.HAPropsSI('Y','P',pressure,'T',temperature_far,'RH',relhum_far)
        mixture_density_far = 1 / HAP.HAPropsSI('Vha','P',pressure,'T',temperature_far,'RH',relhum_far)
        vapor_density_far = mole_fraction_far * mixture_density_far

        enthalpy_vapor_interface = CP.PropsSI('H','T',temperature_interface,'Q',1,'water') \
            - h_ref_water

        # Now determine the film state and transport properties
        enthalpy_film = 0.5 * (enthalpy_interface + enthalpy_far)
        humrat_film = 0.5 * (humrat_interface + humrat_far)
        temperature_film = HAP.HAPropsSI('T','P',pressure,'Hha',enthalpy_film,'HumRat',humrat_film)
        relhum_film = HAP.HAPropsSI('RH','P',pressure,'Hha',enthalpy_film,'HumRat',humrat_film)
        volume = HAP.HAPropsSI('Vha','P',pressure,'Hha',enthalpy_film,'HumRat',humrat_film)
        density = 1 / volume
        c_p = HAP.HAPropsSI('Cha','P',pressure,'Hha',enthalpy_film,'HumRat',humrat_film)
        k = HAP.HAPropsSI('K','P',pressure,'Hha',enthalpy_film,'HumRat',humrat_film)
        mu = HAP.HAPropsSI('mu','P',pressure,'Hha',enthalpy_film,'HumRat',humrat_film)
        nu = mu * volume
        prandtl = c_p * mu / k
        #print(prandtl)
        #print(vapor_density_far, vapor_density_interface, mixture_density_far)

        length_scale = 0.01
        velocity = 1
        nusselt = 1
        h_heat = nusselt * k / length_scale
        stanton_heat = h_heat / (density * velocity * c_p)
        #print(stanton_heat)

        # Mass diffusivity [m2/s]
        # TODO: add a lookup-table with temperature
        diffusivity_mass = 2.5e-5
        schmidt = nu / diffusivity_mass
        lewis = schmidt / prandtl
        #print(lewis)

        stanton_ratio = pow(lewis, 2/3)
        #print(stanton_ratio)

        stanton_mass = stanton_heat / stanton_ratio
        h_mass = stanton_mass * velocity
        #print(h_mass)

        flow_per_area_mass = h_mass * (vapor_density_interface - vapor_density_far)
        flow_per_area_heat = h_heat * (temperature_interface - temperature_far)
        flow_per_area_mass, flow_per_area_heat

        mass_flow_cell = a_cell * flow_per_area_mass
        heat_flow_cell = a_cell * flow_per_area_heat

        # If the state is near saturation and getting closer,
        # then we need to dial down the water vapor influx.
        # Water vapor outflux is okay, as it tends to drive away from saturation.
        if supersaturated_method is 'hack':
            if relhum_film > 1 and mass_flow_cell > 0:
                mass_flow_cell = 0

        water_flow_out = water_flow_in - mass_flow_cell
        humrat_in = humrat_far
        humair_mass_flow_in = dryair_mass_flow_in * (1+humrat_in)
        humair_mass_flow_out = humair_mass_flow_in + mass_flow_cell
        humrat_out = humair_mass_flow_out / dryair_mass_flow_in - 1

        enthalpy_ave_in = enthalpy_far
        total_enthalpy_in = humair_mass_flow_in * enthalpy_ave_in
        total_enthalpy_out = total_enthalpy_in + heat_flow_cell \
            + mass_flow_cell * enthalpy_vapor_interface
        enthalpy_ave_out = total_enthalpy_out / humair_mass_flow_out

        water_total_enthalpy_in = water_flow_in * CP.PropsSI('H','T',temperature_interface,'Q',0,'water')
        water_total_enthalpy_out = water_total_enthalpy_in - heat_flow_cell \
            - mass_flow_cell * enthalpy_vapor_interface
        water_enthalpy_out = water_total_enthalpy_out / water_flow_out
        water_temperature_out = CP.PropsSI('T','H',water_enthalpy_out,'P',pressure,'water')
        #print(water_temperature_out)
        interface_humrat_out = HAP.HAPropsSI('HumRat','P',pressure,'T',water_temperature_out,'RH',relhum_interface)

        t_out = HAP.HAPropsSI('T','P',pressure,'Hha',enthalpy_ave_out,'HumRat',humrat_out)
        rh_out = HAP.HAPropsSI('RH','P',pressure,'Hha',enthalpy_ave_out,'HumRat',humrat_out)

        temperature_far_in = t_out
        relhum_far_in = rh_out
        water_flow_in = water_flow_out
        temperature_interface = water_temperature_out

        nodal_temperature.append(t_out)
        nodal_relhum.append(rh_out)
        nodal_humrat.append(humrat_out)
        nodal_water_flow.append(water_flow_out)
        nodal_water_temperature.append(water_temperature_out)
        nodal_interface_humrat.append(interface_humrat_out)
    
    result = pandas.DataFrame.from_dict(dict(
        T_air=nodal_temperature,
        RelHum=nodal_relhum,
        HumRat_air=nodal_humrat,
        WaterFlow=nodal_water_flow,
        T_water=nodal_water_temperature,
        HumRat_interface=nodal_interface_humrat
        ))
    return result

def plot_results(my_dataframe):
    import matplotlib.pyplot
    import CoolProp.Plots
    import CoolProp.Plots.PsychChart
    import CoolProp.Plots.PsychScript
    #CoolProp.Plots.TwoStage('water',1,Te=280,Tc=300,DTsh=2,DTsc=2,eta_oi=0.9,f_p=0.1,Tsat_ic=290,DTsh_ic=2)
    matplotlib.pyplot.close('all')
    if True:
        from CoolProp.HumidAirProp import HAPropsSI
        from CoolProp.Plots.Plots import InlineLabel
        #from CoolProp.Plots.Common import BasePlot
        #bp = BasePlot()
        #InlineLabel = bp.inline_label

        p = 101325
        Tdb = numpy.linspace(-10,40,100)+273.15

        # Make the figure and the axes
        fig=matplotlib.pyplot.figure(figsize=(10,8))
        ax=fig.add_axes((0.1,0.1,0.85,0.85))

        # Saturation line
        w = [HAPropsSI('W','T',T,'P',p,'R',1.0) for T in Tdb]
        ax.plot(Tdb-273.15,w,lw=2)

        # Humidity lines
        RHValues = [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        for RH in RHValues:
            w = [HAPropsSI('W','T',T,'P',p,'R',RH) for T in Tdb]
            ax.plot(Tdb-273.15,w,'r',lw=1)

        # Humidity lines
        for H in [-20000, -10000, 0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000]:
            #Line goes from saturation to zero humidity ratio for this enthalpy
            T1 = HAPropsSI('T','H',H,'P',p,'R',1.0)-273.15
            T0 = HAPropsSI('T','H',H,'P',p,'R',0.0)-273.15
            w1 = HAPropsSI('W','H',H,'P',p,'R',1.0)
            w0 = HAPropsSI('W','H',H,'P',p,'R',0.0)
            ax.plot(numpy.r_[T1,T0],numpy.r_[w1,w0],'r',lw=1)

        ax.set_xlim(Tdb[0]-273.15,Tdb[-1]-273.15)
        ax.set_ylim(0,0.03)
        ax.set_xlabel(r"$T_{db}$ [$^{\circ}$C]")
        ax.set_ylabel(r"$W$ ($m_{w}/m_{da}$) [-]")

        xv = Tdb #[K]
        for RH in [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
            yv = numpy.array([HAPropsSI('W','T',T,'P',p,'R',RH) for T in Tdb])
            #y = HAPropsSI('W','P',p,'H',65000.000000,'R',RH)
            y = HAPropsSI('W','P',p,'H',45000.000000,'R',RH)
            T_K,w,rot = InlineLabel(xv, yv, y=y, fig=fig, axis = ax)
            string = r'$\phi$='+'{s:0.0f}'.format(s=RH*100)+'%'
            bbox_opts = dict(boxstyle='square,pad=0.0',fc='white',ec='None',alpha = 0.9)
            ax.text(T_K-273.15,w,string,rotation = rot,ha ='center',va='center',bbox=bbox_opts)

    ax.plot(my_dataframe.T_air-273.15,
            my_dataframe.HumRat_air,
            'k.-')
    ax.plot(my_dataframe.T_water-273.15,
            my_dataframe.HumRat_interface,
           'm.-')
    matplotlib.pyplot.show()
    