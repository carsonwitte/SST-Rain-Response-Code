"""
Ocean_Skin_Model.py

Python translation of H. Bellenger's Prognostic Model for SST and SSS including cool-skin, warm-layer, and rain physics. Includes option for Witte modification to remove volume fraction scaling of rain sensible heat flux in cool-skin calculation. 

Call ocean_skin with input vectors of arbitrary length and it will loop through the physics for you. This code could be significantly improved - in speed via vectorization, and in readability via implementation of more specific functions and better variable names.

C.R. Witte
2022
"""

import numpy as np
import pandas as pd
import xarray as xr
import datetime as datetime
import glob
import scipy
import tqdm
import gsw

###########################################################################################

# GLOBAL VARIABLES
cpa = 1004.67  # specific heat of dry air, in J / kg / K (Businger 1982)
cpw = 4000.    # specific heat of liquid water, in J / kg / K (DO BETTER WITH GSW?)
grav = 9.80665 # acceleration due to gravity, m s-2 (given as 9.780326772 by Bellenger, which is the "gamma" parameter in a latitude scaling of g (see COARE function grv))
rhow = 1022.   # density of liquid water, in kg / m3 (DO BETTER WITH GSW?)
rgas = 287.1   # specific ideal gas constant for dry air, in J / kg / K
von = 0.4      # von Karman's "constant"
eps_w = 0.62197 # molecular mass of water over molecular mass of dry air (Gill 1982 k0829, equation 3.1.13)
beta = 0.756 / 1023.343  #Salinity expansion coefficient. Derivative with respect to practical salinity, not mass fraction of salt. Value at 25 Celsius degrees, zero pressure, practical salinity 35, from Gill (1982 k0829, table A.3.1).
depth = 3. # nominal diurnal warm layer and fresh water lens depth, in m (Zeng and Beljaars 2005)


###########################################################################################

def esat(t, p):
    '''
    Computes saturation vapor pressure of water in Pa. From Buck, 1981, J. Appl. Meteor. 20, 1527-1532, equation (8).
    
    Inputs:
        t    - temperature [K]
        p    - air pressure [Pa]
    '''
    esat = (1.0007 + 3.46e-8 * p) * 6.1121e2 * np.exp(17.502 * (t - 273.15) / (t - 32.18))
    return esat

###########################################################################################

def therm_expans(t):
    '''
    Calculate alpha, thermal expansion coefficient for seawater. Could be replaced with the GSW toolbox function
    
    Input:
        t    - temperature [K]
    '''
    alpha = 2.1e-5 * (t - 269.95)**0.79
    return alpha

###########################################################################################

def phiw(zL):
    '''
    Calculate nondimensional temperature gradient function (wrote a better vectorized version of this below)
    
    Input:
        zL  - z over L
    '''
    if zL < 0:
        phiw = (1. - 16. * zL)**(- 0.5)
    else:
        phiw = 1. + (5. * zL + 4. * zL**2)/(1. + 3. * zL + 0.25 * zL**2) #(Takaya et al. 2010)
    
    return phiw

###########################################################################################

def fV(z, rain):
    '''
    Calculate fraction of rain volume entering the ocean and deposited within depth "z"
    Adapted from H. Bellenger 2016 Fortran code
    
    Inputs:
        z     - depth [m] (z < 0)
        rain  - rain mass flux [kg/m^2/s]
    '''
    a = 100 #coefficient of proportionality between raindrop radius and penetration depth (Manton, 1973)

    # Schlussel et al. 1997, Table 1:
    ζ = np.array([0.0026091, 0.0022743, 0.0015406, 0.0012281, 0.0008795, 0.00077123, 0.00057451, 0.000438, 6.7228e-5, 6.4955e-5, 4.4234e-5, 3.3906e-5, 2.7433e-6, 4.0283e-7])
    ψ = np.array([0.60107, 0.29968, 0.5563, 1.80858, 0.2175, 0.33961, 0.96368, 0.65081, 0.5967, 2.7661, 2.2812, 2.7674, 2.6095, 6.5308])

    if rain > 0:
        rc = 0.4 #[mm]
        z_mm = z * 1000. #[mm] and <0
        Λ = 4.1 * (rain * 3600.)**(- 0.21) # mm-1
        f0 = (1 + 2*Λ*rc + 0.5*(2*Λ*rc)**2 + (1/6)*(2*Λ*rc)**3) * np.exp(-2*Λ*rc) 
        dfv = np.sum(ζ * np.exp(-ψ*Λ*abs(z_mm)/100))
        fV = -dfv*abs(z_mm) + f0
    else:
        fV = 0
    
    return fV

###########################################################################################

def sens_heat_rain(rain, t, q, rhoa, xlv, t_int, p):
    '''
    Computes heat flux due to rainfall, in W m-2, positive upward. From Gosnell et al., 1995 equations 11 and 12
    Modified by C. Witte from original Bellenger formulation, which used interface temperature/humidity instead of air temperature/humidity in the Clausius-Clapeyron step
    
    Inputs:
        rain      - rainfall in kg m-2 s-1
        t         - air temperature, in K
        q         - specific humidity, in kg/kg
        rhoa      - density of moist air  (kg / m3)
        xlv       - latent heat of evaporation (J / kg)
        t_int     - air-sea interface temperature, in K
        p         - surface pressure, in Pa
    '''

    es = esat(t_int, p) # saturation pressure of wator vapor, in Pa, reduced for salinity, Kraus 1972 page 46
    q_int = eps_w * (es / (p - (1. - eps_w) * es)) #specific (saturation) humidity at ocean interface
    wetc = eps_w * xlv * q / (rgas * t**2) #this is the part that should be using air temperature & humidity
    dwat = 2.11e-5 * (t / 273.15)**1.94 #water vapour diffusivity 
    t_celsius = t - 273.15 
    dtmp = (1. + 3.309e-3 * t_celsius - 1.44e-6 * t_celsius**2) * 0.02411 / (rhoa * cpa) #heat diffusivity 

    # Gosnell 1995 k0991, equation (11):
    alfac =  1. / (1. + (wetc * xlv * dwat) / (cpa * dtmp)) # wet bulb factor 

    # Gosnell 1995 k0991, equation (12):
    sens_heat_rain =  rain * alfac * cpw * (t_int - t + (q_int - q) * xlv / cpa)
    
    return sens_heat_rain

###########################################################################################

def sens_heat_rain_orig(rain, t, q, rhoa, xlv, t_int, p):
    '''
    ORIGINAL VERSION FROM BELLENGER PORT BEFORE MODIFICATION
    Computes heat flux due to rainfall, in W m-2, positive upward. From Gosnell et al., 1995 equations 11 and 12
    
    Inputs:
        rain      - rainfall in kg m-2 s-1
        t         - air temperature, in K
        q         - specific humidity, in kg/kg
        rhoa      - density of moist air  (kg / m3)
        xlv       - latent heat of evaporation (J / kg)
        t_int     - air-sea interface temperature, in K
        p         - surface pressure, in Pa
    '''

    es = esat(t_int, p) * 0.98 # saturation pressure of wator vapor, in Pa, reduced for salinity, Kraus 1972 page 46
    q_int = eps_w * (es / (p - (1. - eps_w) * es)) #specific (saturation) humidity at ocean interface
    wetc = eps_w * xlv * q_int / (rgas * t_int**2)
    dwat = 2.11e-5 * (t / 273.15)**1.94 #water vapour diffusivity IDENTICAL
    t_celsius = t - 273.15 
    dtmp = (1. + 3.309e-3 * t_celsius - 1.44e-6 * t_celsius**2) * 0.02411 / (rhoa * cpa) #heat diffusivity IDENTICAL

    # Gosnell 1995 k0991, equation (11):
    alfac =  1. / (1. + (wetc * xlv * dwat) / (cpa * dtmp)) # wet bulb factor DIFFERENT - uses wetc instead of dqs_dt

    # Gosnell 1995 k0991, equation (12):
    sens_heat_rain =  rain * alfac * cpw * (t_int - t + (q_int - q) * xlv / cpa)
    
    return sens_heat_rain

###########################################################################################

def mom_flux_rain(u, rain):
    '''
    Compute momentum flux due to rainfall [Pa]
    
    Inputs:
        u      - difference of velocity between air and sea, including gustiness [m/s]
        rain   - rain mass flux [kg/m^2/s]
    '''
    
    mom_flux_rain = 0.85 * u * rain   # (Caldwell 1971 k1001, equation (1) and 15 % reduction in speed as in paragraph 3.a, maybe adequate if u is the wind at 10 m)
    
    return mom_flux_rain

###########################################################################################

def microlayer(tkt, tks, hlb, tau, s_subskin, al, xlv, taur, rf, rain, qcol, fV_flag):
    '''
    Compute skin temperature and salinity given subskin values and approximate thicknesses of the layers; return skin values and revised thicknesses
    
    Inputs:
        tkt       - approximate thickness of thermal microlayer (cool skin) [m]
        tks       - approximate thickness of haline microlayer [m] 
        hlb       - latent heat flux at the surface, positive upward [W/m^2]
        tau       - wind stress, turbulent part only [Pa]
        s_subskin - subskin salinity [ppt]
        al        - thermal expansion coefficient [1/K]
        xlv       - latent heat of evaporation [J/kg]
        taur      - momentum flux due to rainfall [Pa]
        rf        - sensible heat flux at the surface due to rainfall, positive upward [W/m^2]
        rain      - rain mass flux [kg/m^2/s]
        qcol      - net flux at the surface (without sensible heat flux due to rain) [W/m^2]
        fV_flag   - if True, scale rain sensible heat flux by volume fraction as specified in Bellenger 2017 when calculating skin temp
    
    Outputs:
        dter      - delta temperature in diffusive thermal microlayer (air-sea interface minus subskin temperature) [C] (aka K cause it's a difference)
        dser      - delta salinity in diffusive haline microlayer (air-sea interface minus subskin salinity) [ppt]
        tkt       - revised thickness of thermal microlayer (cool skin) [m]
        tks       - revised thickness of haline microlayer [m] 
    '''

    #real, dimension(size(qcol)):: usrk, usrct, usrcs, alq
    #real xlamx(size(qcol)) ! Saunders coefficient
    
    visw = 1e-6  # viscosity of water?
    tcw = 0.6    # thermal conductivity of water
    mu = 0.0129e-7 # molecular salinity diffusivity [m^2/s], Kraus and Businger, page 47
    kappa = 1.49e-7 # thermal diffusivity [m^2/s]

    # a and b coefficients for the power function fitting the TKE flux carried by rain: Fk = a * R**b, derived from the exact solution of Soloviev and Lukas 2006 (Schlussel et al 1997, Craeye and Schlussel 1998)
    afk = 4e-4
    bfk = 1.3

    usrk = (afk / rhow)**(1/3) * (rain * 3600)**(bfk / 3) #Equivalent friction velocity due to the TKE input by the penetrating raindrops Fk
    
    #Buoyancy Flux with option whether to scale rain sensible heat flux  & precip mass flux by the volume fraction of rainwater within microlayer
    if fV_flag:
        alq = al * (qcol + rf * (1 - fV(tkt, rain))) - beta * s_subskin * cpw * (hlb / xlv - rain * (1 - fV(tks, rain)))
    else:
        alq = al * (qcol + rf) - beta * s_subskin * cpw * (hlb / xlv - rain)

    # Friction velocities in the air:
    usrct = np.sqrt((tau + (1. - fV(tkt, rain)) * taur) / rhow + (fV(0., rain) - fV(tkt, rain)) * usrk**2)
    usrcs = np.sqrt((tau + (1. - fV(tks, rain)) * taur) / rhow + (fV(0., rain) - fV(tks, rain)) * usrk**2)
  

    if alq > 0:
        #Fairall 1996 982, equation (14):
        xlamx = 6. * (1. + (16. * grav * cpw * rhow * visw**3 * alq / (tcw**2 * usrct**4 ))**0.75)**(- 1. / 3.)

        #Fairall 1996 982, equation (12):
        tkt = xlamx * visw / usrct
        
        # From Saunders 1967 (4)
        tks = xlamx * mu * (kappa / mu)**(2. / 3.) * visw * cpw * rhow / (tcw * usrcs)
       
    else:
        xlamx = 6 # prevent excessive warm skins
        tkt = min(.01, xlamx * visw / usrct) # Limit tkt
        tks = min(.001, xlamx * mu * (kappa / mu)**(2. / 3.) * visw * cpw * rhow / (tcw * usrcs))

    # Fairall 1996 982, equation (13) with option whether to scale rain sensible heat flux by the volume fraction of rainwater within microlayer
    if fV_flag:
        dter = - (qcol + rf * (1 - fV(tkt, rain))) * tkt / tcw
    else:
        dter = - (qcol + rf) * tkt / tcw

    dser = s_subskin * (hlb / xlv - rain * (1 - fV(tks, rain))) * tks / (rhow * mu) # eq. fresh skin

    return dter, dser, tkt, tks

###########################################################################################

def near_surface(ds_ns, dt_ns, tau, taur, hlb, rhoa, xlv, dtime, t_ocean_1, s1, rain, q_pwp, depth_1, spread_flag):
    '''
    Calculate heat and freshwater budget for near-surface layer of the ocean (get subskin T & S from foundation T & S plus fluxes)
    
    Inputs:
        tau        - turbulent wind stress at the surface [Pa]
        taur       - momentum flux due to rainfall [Pa]
        hlb        - turbulent latent heat flux [W/m^2]
        rhoa       - density of moist air  [kg/m^3]
        xlv        - latent heat of evaporation [J/kg]
        dtime      - time step [s]
        t_ocean_1  - input sea temperature, at depth_1 [K]
        s1         - salinity at depth_1 [ppt]
        rain       - rain mass flux [kg/m^2/s]
        q_pwp      - net flux absorbed by the warm layer (part of the solar flux absorbed at "depth"), minus surface fluxes, in [W m^2]
        depth_1    - depth of t_ocean_1 and s1 measurements
        spread_flag - if True, include spreading term to account for temporal evolution of the layer
        
    Outputs:
        al         - thermal expansion coefficient [1/K]
        ds_ns      - "delta salinity near surface" [ppt]. Salinity variation in the near-surface turbulent layer. That is subskin salinity minus foundation salinity.
        dt_ns      - "delta temperature near surface" [K aka C]. Temperature variation in the near-surface turbulent layer. That is subskin temperature minus foundation temperature.
        t_subskin  - subskin temperature, in K
        s_subskin  - subskin salinity, in ppt
    '''
    khor = 1. / 1.5e4  # Parameter for the lens spread, in m-1. Inverse of the size of the lens.
    umax = 15.
    fact = 1.

    # Temperature and salinity profiles change with wind:
    u = 28. * np.sqrt(tau / rhoa)

    # assign value for eta based on wind speed and sign of temperature gradient 
    # just using if statements for now as we're naively looping line by line
    if dt_ns < 0:
        if u >= umax:
            eta = 1. / fact
        elif u <= 2:
            eta = 2. / (fact * umax)
        else: #for 2 < u < umax
            eta = u / (fact * umax)
    else: #outside of fresh layers, leave eta constant
        eta = 0.3
        
    if depth_1 < depth: #if the depth of the measurements is within the nominal warm/fresh layer thickness, scale the "foundation" temperature to the assumed depth
        correction = 1. - (depth_1 / depth)**eta #(neglecting microlayer thickness compared to depth_1 and depth)
        t_fnd = t_ocean_1 - dt_ns * correction
        s_fnd = s1 - ds_ns * correction
    else: #otherwise assume well-mixed below the foundation depth
        t_fnd = t_ocean_1
        s_fnd = s1

    al = therm_expans(t_fnd)

    # Bellenger 2017 k0976, equation (13):
    buoyf = al * grav / (rhow * cpw) * q_pwp - beta * s_fnd * grav * (hlb / xlv - rain) / rhow

    usrc = np.sqrt((tau + taur) / rhow)
    drho = rhow * (- al * dt_ns + beta * ds_ns)

    # Case of stable stratification and negative flux, Bellenger 2017  k0976, equation (15):
    if buoyf < 0 and drho < 0:
        buoyf = np.sqrt(- eta * grav / (5. * depth * rhow) * drho) * usrc**2
    elif buoyf == 0:
        buoyf = np.finfo(np.float64).tiny #smallest possible non-zero positive number to avoid divide-by-zero issues

    Lmo = usrc**3 / (von * buoyf)

    # Equation (14) for temperature. Implicit scheme for time integration:
    # \Delta T_{i + 1} - \Delta T_i = \delta t (Bt + At \Delta T_{i + 1})
    At = - (eta + 1.) * von * usrc / (depth * phiw(depth / Lmo))

    # Lens horizontal spreading:
    if drho < 0 and ds_ns < 0 and spread_flag:
        At = At - (eta + 1.) * khor * np.sqrt(depth * grav * abs(drho) / rhow)

    Bt = q_pwp / (depth * rhow * cpw * eta / (eta + 1.))
    dt_ns = (dtime * Bt + dt_ns) / (1 - dtime * At)
    
    #try an alternative implicit integration: 
    # \Delta T_{i + 1} - \Delta T_i = \delta t (Bt + At \Delta T_{i})
    #dt_ns = (Bt + At*dt_ns)*dtime + dt_ns

    # Equation (14) for salinity:
    # \frac{\partial \Delta S}{\partial t} = (\Delta S + S_\mathrm{fnd}) B_S + A_S \Delta S
    As = - (eta + 1.) * von * usrc / (depth * phiw(depth / Lmo))

    # Lens horizontal spreading:
    if drho < 0 and ds_ns < 0 and spread_flag:
        As = As - (eta + 1.) * khor * np.sqrt(depth * grav * abs(drho) / rhow)

    Bs = (hlb / xlv - rain) * (eta + 1.) / (depth * rhow * eta)

    # Implicit scheme for time integration:
    ds_ns = (dtime * Bs * s_fnd + ds_ns) / (1 - dtime * (As + Bs))

    t_subskin = t_fnd + dt_ns
    s_subskin = s_fnd + ds_ns
    
    return al, dt_ns, ds_ns, t_subskin, s_subskin

###########################################################################################

def bulk_flux(u, t_ocean_1, s1, rain, hf, hlb, rnl, tau, rhoa, xlv, rf, dtime, rns, ds_ns, dt_ns, depth_1, jcool, jwarm, rain_effect, fV_flag, spread_flag):
    '''
    Bulk flux loop similar to COARE3.5 but with rain physics added
    
    Inputs:
        u           - wind speed relative to the sea surface [m/s] (i. e. taking current vector into account)
        t_ocean_1   - input sea temperature, at depth_1 [K]
        s1          - salinity at depth_1 [ppt]
        rain        - rain mass flux [kg/m^2/s]
        hf          - turbulent sensible heat flux, positive upward [W/m^2]
        hlb         - latent heat flux at the surface, positive upward [W/m^2]
        rnl         - net longwave radiation, positive upward [W/m^2]
        tau         - wind stress, turbulent part only [Pa]
        rhoa        - density of moist air  [kg/m^3]
        xlv         - latent heat of evaporation (J / kg)
        rf          - sensible heat flux at the surface due to rainfall, positive upward [W/m^2]
        dtime       - time step [s]
        rns         - net shortwave radiation [W/m^2]
        dt_ns       - delta temperature near surface (subskin minus foundation temperature) [C] (aka K)
        ds_ns       - delta salinity near surface (subskin minus foundation salinity) [ppt]
        depth_1     - depth of bulk measurement [m]
        jcool       - boolean flag to turn on skin-layer physics
        jwarm       - boolean flag to turn on warm-layer physics
        rain_effect - boolean flag to turn on rain physics
        fV_flag     - if True, scale rain sensible heat flux by volume fraction as specified in Bellenger 2017 when calculating skin temp
        spread_flag - if True, include spreading term to account for temporal evolution of the layer

    Outputs:
        tkt       - thickness of thermal microlayer (cool skin) [m]
        tks       - thickness of haline microlayer [m] 
        taur      - momentum flux due to rain [Pa]
        dter      - delta temperature in diffusive thermal microlayer (air-sea interface minus subskin temperature) [C] (aka K cause it's a difference)
        dser      - delta salinity in diffusive haline microlayer (air-sea interface minus subskin salinity) [ppt]
        t_int     - interface temperature, [K]
        s_int     - interface salinity [ppt]
        dt_ns     - updated delta temperature near surface (subskin minus foundation temperature) [C] (aka K cause it's a difference)
        ds_ns     - updated delta salinity near surface (subskin minus foundation salinity) [ppt]

    '''
    
    #Soloviev solar absorption profile - original from Bellenger code:
    #fxp = 1. - (0.28 * 0.014 + 0.27 * 0.357 * (1. - np.exp(- depth / 0.357)) + 0.45 * 12.82 * (1.- np.exp(- depth / 12.82))) / depth
    #modified to look the same as COARE warm-layer code - was missing the first exponent term...
    fxp = 1. - (0.28*0.014*(1-np.exp(- depth / 0.014)) + 0.27*0.357*(1-np.exp(- depth / 0.357)) + 0.45*12.82*(1.- np.exp(- depth / 12.82))) / depth
    
    #Soloviev solar absorption profile - taken by C.Witte from Kudryavtsev & Soloviev 1990 pg. 623 Eq 1 - not as good
    #fxp = (0.28*np.exp(-71.5*depth) + 0.27*np.exp(2.8*depth) + 0.45*np.exp(0.07*depth))

    
    tau_0 = 1e-3 # in N m-2

    #calculate rain momentum flux if rain physics are turned on, or set to zero if they're not
    if rain_effect:
        taur = mom_flux_rain(u, rain)
    else: #Bellenger had another piece of this if statement just setting up a zero value... I don't see the use: "if (jwarm .or. jcool) null_array = 0."
        taur = 0.

    tau_with_min = tau + tau_0 * (1 - np.exp(- tau_0 / tau))
    
    #-----NEAR SURFACE LAYER-----
    if jwarm:
        if rain_effect:
            #q_pwp = fxp * rns - (hf + hlb + rnl + rf)
            q_pwp = fxp * rns - (hf + hlb + rnl + rf)
            [al, dt_ns, ds_ns, t_subskin, s_subskin] = near_surface(ds_ns, dt_ns, tau, taur, hlb, rhoa, xlv, dtime, t_ocean_1, s1, rain, q_pwp, depth_1, spread_flag)
        else:
            #q_pwp = fxp * rns - (hf + hlb + rnl)
            #solar flux is already positive upward, so try switching its sign...
            q_pwp = fxp * rns - (hf + hlb + rnl)
            [al, dt_ns, ds_ns, t_subskin, s_subskin] = near_surface(ds_ns, dt_ns, tau, taur, hlb, rhoa, xlv, dtime, t_ocean_1, s1, 0, q_pwp, depth_1, spread_flag)
    else:
        al = therm_expans(t_ocean_1)
        t_subskin = t_ocean_1
        s_subskin = s1
        
    #-----SKIN LAYER-----
    if jcool:
        # First guess:
        tkt = 0.001
        tks = 5e-4

        for i in np.arange(0,3): #just calls the microlayer loop three times, updating the input parameters based on the outputs (usually converges after 1 iteration)
            # Cool skin
            dels = rns * (0.065 + 11. * tkt - 6.6e-5 / tkt * (1. - np.exp(- tkt / 8e-4))) # equation 16 Ohlmann
            qcol = rnl + hf + hlb - dels
            if rain_effect:
                [dter, dser, tkt, tks] = microlayer(tkt, tks, hlb, tau_with_min, s_subskin, al, xlv, taur, rf, rain, qcol, fV_flag)
            else:
                [dter, dser, tkt, tks] = microlayer(tkt, tks, hlb, tau_with_min, s_subskin, al, xlv, taur, 0, 0, qcol, fV_flag)
    else:
        tkt = 0.
        tks = 0.
        dter = 0.
        dser = 0.

    t_int = t_subskin + dter
    s_int = s_subskin + dser
    
    return tkt, tks, taur, dter, dser, t_int, s_int, ds_ns, dt_ns

###########################################################################################

def ocean_skin(time, u, t_ocean_1, s1, t, q, rsds, p, rain, hf, hlb, rnl, tau, depth_1, jcool, jwarm, rain_effect, fV_flag, spread_flag=True):
    '''
    Calculate ocean skin temperature including cool-skin, warm layer, and rain physics following Bellenger et al., 2016
    
    Inputs:
        time      - elapsed time in seconds
        u         - wind speed relative to the sea surface, i. e. taking current vector into account. In m s-1.
        t_ocean_1 - input sea temperature, at depth_1, in C
        s1        - salinity at depth_1, in ppt
        t         - air temperature, in C
        q         - specific humidity, in kg/kg
        rsds      - surface downwelling shortwave radiation, positive upward, in W / m2
        p         - surface pressure, in hPa (aka mbar)
        rain      - rainfall in mm/hr
        hf        - turbulent sensible heat flux, in W m-2
        hlb       - latent heat flux at the surface (W m-2)
        rnl       - net longwave radiation, positive upward, in W m-2
        tau       - wind stress, turbulent part only, in Pa
        depth_1   - depth of bulk measurement
        jcool       - boolean flag to turn on skin-layer physics
        jwarm       - boolean flag to turn on warm-layer physics
        rain_effect - boolean flag to turn on rain physics
        fV_flag     - if True, scale rain sensible heat flux by volume fraction as specified in Bellenger 2017 when calculating skin temp
        spread_flag - if True, include spreading term to account for temporal evolution of the layer
   
    Outputs:
        t_int     - air-sea interface temperature, in C
        s_int     - air-sea interface salinity, in ppt
        tkt       - thermal molecular diffusion microlayer (cool skin) thickness in m
        tks       - haline molecular diffusion microlayer thickness in m
        dter      - delta temperature in diffusive thermal microlayer (air-sea interface minus subskin temperature) in C
        dser      - delta salinity in diffusive haline microlayer (air-sea interface minus subskin salinity) in ppt
        dt_ns     - delta temperature near surface (subskin minus foundation temperature) in C
        ds_ns     - delta salinity near surface (subskin minus foundation salinity) in ppt
        rf        - rain sensible heat flux, in W m-2
        taur      - momentum flux due to rain, in Pa

    '''
    #----INITIAL VALUES----
    t_ocean_1 = t_ocean_1 + 273.15 # Celsius degrees to K
    t_int = t_ocean_1[0]              # initial interface temp guess by setting to input temp
    ds_ns = 0.
    dt_ns = 0.
    old_time = 0.  # in s
    
    #----CONVERSION AND SETUP---
    albedo = 0.055 # daily average
    
    t = t + 273.15 # Celsius degrees to K
    p = p * 100.   # convert hPa to Pa
    rain  = rain / 3600. # convert mm/h -> kg m-2 s-1
    rhoa = p / (rgas * t * (1. + (1. / eps_w - 1.) * q)) #calculate density of moist air (kg / m3)

    #the following is structured as a loop in fortran that goes line by line of data (long-term I would prefer to pass full arrays rather than looping)
    #initialize a bunch of empty lists to fill as we go through the loop (terrible speed-wise but I'm trying to reproduce as exactly as possible before optimizing)
    tkt_out = []
    tks_out = []
    taur_out = []
    dter_out = []
    dser_out = []
    t_int_out = []
    s_int_out = []
    ds_ns_out = []
    dt_ns_out = []
    rf_out = []
    
    for idx in tqdm.tqdm(np.arange(len(time))):
        xlv = (2.501 - 0.00237 * (t_int - 273.15)) * 1e6  #calculate latent heat of evaporation for this t_int (J / kg)
        
        rf = sens_heat_rain(rain[idx], t[idx], q[idx], rhoa[idx], xlv, t_int, p[idx]) #Gosnell 95 rain sensible heat flux

        dtime = time[idx] - old_time
        rns = (1. - albedo) * rsds[idx]
        
        #note that ds_ns and dt_ns are both inputs and outputs to this function, and will thus update every loop
        [tkt, tks, taur, dter, dser, t_int, s_int, ds_ns, dt_ns] = bulk_flux(u[idx], t_ocean_1[idx], s1[idx], rain[idx], hf[idx], hlb[idx], rnl[idx], tau[idx], rhoa[idx], xlv, rf, dtime, rns, ds_ns, dt_ns, depth_1, jcool, jwarm, rain_effect, fV_flag, spread_flag)

        tkt_out.append(tkt)
        tks_out.append(tks)
        taur_out.append(taur)
        dter_out.append(dter)
        dser_out.append(dser)
        t_int_out.append(t_int)
        s_int_out.append(s_int)
        ds_ns_out.append(ds_ns)
        dt_ns_out.append(dt_ns)
        rf_out.append(rf)

        old_time = time[idx] #update the "old_time" for the next iteration
    
    return np.array(t_int_out) - 273.15, np.array(s_int_out), np.array(tkt_out), np.array(tks_out), np.array(dter_out), np.array(dser_out), np.array(dt_ns_out), np.array(ds_ns_out), np.array(rf_out), np.array(taur_out)
