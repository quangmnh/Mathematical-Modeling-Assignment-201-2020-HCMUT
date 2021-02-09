''' Dependencies '''
import numpy as np
import math

''' Supported functions'''

def convert_ppm_to_mgm3(ppm):
    return ppm/0.554

def convert_mgm3_to_ppm(mgm3):
    return mgm3*0.554

''' Define constants '''
eta_heatCO2 = 0.057
A_flr = 1.4*10**4
A_roof = 0.1 * A_flr
A_side = 0 * A_flr
P_blow = 0
phi_pad = 0
g = 9.81
phi_extCO2 = 7.2*10**4
K_thscr = 0.05*10**-3
M_air = 28.96
h_elevation = 0
h_elevation_top = 3.8
h_elevation_mean_air = 4.2
cap_CO2_air = h_elevation_top
cap_CO2_top = h_elevation_mean_air - h_elevation_top
rho_air0 = 1.2
C_d_Gh = 0.75
C_w_Gh = 0.09
zeta_insscr = 1
c_leakage = 10**-4
eta_roof_thr = 0.9
M_CH2O = 30*10**-3
C_max_buf = 20*10**-3
R = 8.314*10**3
S = 710
eta_pad = 0
LAI = 2.5
J_MAX_25leaf = 210
eta_shscrCd = 0
eta_shscrCw = 0
c_Gamma = 1.7
eta_CO2_air_stom = 0.67
Ej = 37*10**3
T_25 = 298.15
H = 22*10**4
alpha = 0.385
PAR_can = 100
Theta = 0.7
h_side_roof = 0
eta_side_thr = eta_roof_thr
h_roof = 0.68
phi_vent_forced = 0
CO2_out = 668

T_air_out = 1       # T_air - T_out = 1
T_top_out = 1       # T_top - T_out = 1

''' Define variables which value is depended on time'''
U_blow = 0.2
U_extCO2 = 0.1
U_thscr = 0.5
U_vent_forced = 0.5
U_pad = 0.5
U_roof = abs(28.4 - 1.6) / 100
U_side = 0.5
U_shscr = 0.5


v_wind = 6.3
T_out = 17.8
T_air = 21.8999999966472
T_top = T_air + 1
T_can = T_air + 15
C_buf = 20*10**-3

def MC_blow_air(eta_heatCO2, U_blow, P_blow, A_flr):
    '''
    Calculate the flow of CO2 from blower to the greenhouse air

    :param eta_heatCO2: Amount of CO2 generated when 1J of sensible heat is generated
    by the direct air heater
    :param U_blow: The control valve of the direct air heater
    :param P_blow: The heat capacity of the direct air heater
    :param A_flr: Greenhouse floor area
    '''
    return (eta_heatCO2*U_blow*P_blow)/A_flr
#------------------------------------------------------------------


def MC_ext_air(U_extCO2, phi_extCO2, A_flr):
    '''
    Calculate the CO2 supply rate (from the provider)

    :param U_extCO2: The control valve of the external CO2 source
    :param phi_extCO2: The capacity of the external CO2 source
    :param A_flr: Greenhouse floor area
    '''
    return (U_extCO2*phi_extCO2)/A_flr
#------------------------------------------------------------------


def MC_pad_air(U_pad, phi_pad, A_flr, CO2_out, CO2_air):
    '''
    Calculate the flow of CO2 from the pad and fan system to greenhouse air

    :param U_pad: Pad and fan control
    :param phi_pad: Capacity of the air flux through the pad and fan system
    :param CO2_out: Outdoor CO2 concentration
    :param CO2_air: CO2-concentration of the greenhouse air
    '''
    return ((U_pad*phi_pad)/A_flr)*(CO2_out - CO2_air)
#------------------------------------------------------------------


def MC_air_top(f_thscr, CO2_air, CO2_top):
    '''
    Calculate the flow of CO2 from the greenhouse air to the top compartment air

    :param f_thscr: Air flux rate through the thermal screen
    :param CO2_air: CO2-concentration of the greenhouse air
    :param CO2_top: CO2-concentration of the top compartment air
    '''
    return f_thscr*(CO2_air - CO2_top)

def f_thscr(U_thscr, K_thscr, T_air, T_top, g, rho_mean_air, rho_air, rho_top):
    '''
    Calculate the air flux rate through the thermal screen

    :param U_thscr: Control of the thermal screen
    :param K_thscr: The thermal screen flux coefficient
    :param T_air: Greenhouse air temperature
    :param T_top: Greenhouse top compartment air temperature
    :param g: Gravitational acceleration
    :param rho_mean_air: The mean density of the greenhouse and the top compartment air
    :param rho_air: The density of greenhouse air
    :param rho_top: The density of top compartment air
    '''
    Thscr = U_thscr*K_thscr*(abs(T_air - T_top))**(2/3)
    NotThscr = (1 - U_thscr)*( (g*(1-U_thscr))/(2*rho_mean_air)*abs(rho_air - rho_top) )**(1/2)

    return Thscr + NotThscr

def rho_air(rho_air0, g, M_air, h_elevation, T_air, R):
    '''
    Calculate the density of the greenhouse air

    :param rho_air0: The density of the air at sea level
    :param g: Gravitational acceleration
    :param M_air: The molar mass of air
    :param h_elevation: The altitude of the greenhouse above sea level
    :param T_air: Greenhouse air temperature
    :param R: The molar gas constant
    '''
    # return rho_air0*math.exp( (g*M_air*h_elevation)/((273.15+T_air)*R) )
    T_air += 273.15
    return 101325*((1-2.25577*10**-5 * h_elevation)**5.25588)/(R*T_air)

def rho_top(rho_air0, g, M_air, h_elevation_top, T_top, R):
    '''
    Calculate the density of the top compartment air

    :param rho_air0: The density of the air at sea level
    :param g: Gravitational acceleration
    :param M_air: The molar mass of air
    :param h_elevation: The altitude of the greenhouse above sea level
    :param T_top: Top compartment temperature
    :param R: The molar gas constant
    '''

    # return rho_air0*math.exp( (g*M_air*h_elevation_top)/((273.15+T_top)*R) )
    T_top += 273.15
    return 101325 * ((1 - 2.25577 * 10 ** -5 * h_elevation_top) ** 5.25588) / (R * T_top)

def rho_mean_air(rho_air0, g, M_air, h_elevation_mean_air, T_air, T_top, R):
    '''
    Calculate the mean density of the greenhouse air and outdoor air

    :param rho_air: The density of the greenhouse air
    :param rho_out: The density of the outdoor air
    '''
    # return rho_air0*math.exp( (g*M_air*h_elevation_mean_air)/((273.15 + (T_air + T_top)/2 )*R) )
    T_air += 273.15
    T_top += 273.25
    return 101325 * ((1 - 2.25577 * 10 ** -5 * h_elevation_mean_air) ** 5.25588) / (R * (T_air+T_top)/2)

#------------------------------------------------------------------


def MC_air_out(f_vent_side, f_vent_forced, CO2_air, CO2_out):
    '''
    Calculate the flow of CO2 from the greenhouse air to the outdoor

    :param f_vent_side: Natural ventilation rate for side window
    :param f_vent_fourced: Forced ventilation rate for side window
    :param CO2_air: CO2-concentration rate of the greenhouse air
    :param CO2_out: CO2-concentration rate of the outdoor
    '''
    return (f_vent_side + f_vent_forced)*(CO2_air - CO2_out)

def f_vent_roof_side(C_d, A_flr, U_roof, U_side, A_roof, A_side, g, h_side_roof,
        T_air, T_out, T_mean_air, C_w, v_wind):
    '''
    Calculate the ventilation rate through both the roof and the side vents

    :param C_d: Discharge coefficient
    :param A_flr: Greenhouse floor area
    :param U_roof: Control of the aperture of the roof vent
    :param U_side: Control of the side ventilators
    :param A_roof: Maximum roof ventilation area
    :param A_side: Maximum sidewall ventilation area
    :param g: Gravitational acceleration
    :param h_side_roof: Vertical distance between mid-points of side wall and
    roof ventilation openings
    :param T_air: Greenhouse air temperature
    :param T_out: Outdoor temperature
    :param T_mean_air: Mean temperature of the greenhouse air and the outside air
    :param C_w: Global wind pressure coefficent
    :param v_wind: The wind speed
    '''
    if (U_side*A_side == 0 and U_roof*A_roof == 0):
        return 0
    else:
        return (C_d/A_flr)* \
            (((U_roof**2)*(U_side**2)*(A_roof**2)*(A_side**2)) / \
             ((U_roof**2)*(A_roof**2)+(U_side**2)*(A_side**2)) \
             * (2*g*h_side_roof*(T_air - T_out)) / (T_mean_air + 273.15) \
             + (((U_roof*A_roof + U_side*A_side)/2)**2) * C_w * v_wind ** 2)**(1/2)

def T_mean_air(T_air, T_out):
    '''
        Calculate the mean temperature

        :param T_air: The temperature of the greenhouse air
        :param T_out: The temperature of the outdoor air
    '''
    return (T_air + T_out) / 2

def C_d(C_d_Gh, eta_shscrCd, U_shscr):
    '''
    Calculate the discharge coefficent

    :param C_d_Gh: the discharge coefficient determined for a greenhouse
    without an external shading screen
    :param eta_shscrCd: parameter that determines the effect of the
    shading screen on the discharge coefficient
    :param U_shscr: Control of the external shading screen
    '''
    return C_d_Gh*(1 - eta_shscrCd*U_shscr)

def C_w(C_w_Gh, eta_shscrCw, U_shscr):
    '''
    Calculate the global wind pressure coefficent

    :param C_w_Gh: the global wind pressure coefficient for a greenhouse
    without an external shading screen
    :param eta_shscrCw: parameter that determines the effect of the
    shading screen on the global wind pressure coefficient
    :param U_shscr: Control of the external shading screen
    '''
    return C_w_Gh*(1 - eta_shscrCw*U_shscr)

def eta_insscr(zeta_insscr):
    '''
    Calculate the ventilation rate's reduce factor created by insect screens

    :param zeta_insscr: The screen porosity
    '''
    return zeta_insscr*(2 - zeta_insscr)

def f_leakage(c_leakage, v_wind):
    '''
    Calculate the leakage rate

    :param c_leakage: The leakage coefficent
    :param v_wind: The wind speed
    '''
    if v_wind < 0.25:
        return 0.25*c_leakage
    else:
        return v_wind*c_leakage

def f_vent_side(eta_insscr, f2_vent_side, f_leakage, U_thscr, f_vent_roof_side,
        eta_side, eta_roof, eta_roof_thr):
    '''
    Calculate the ventilation rate through the side vents

    :param eta_insscr: The ventilation rate's reduce factor created by insect screens
    :param f2_vent_side: The ventilation rate for sidewall ventilation only
    :param f_leakage: The leakage rate
    :param U_thscr: Control of the thermal screen
    :param f_vent_roof_side: The ventilation rate through both the roof and the side vents
    :param eta_side: The ratio between the side vents area and total ventilation area
    :param eta_roof: The ratio between the roof vents area and total ventilation area
    :param eta_roof_thr: The threshold value for which there is no chimney effect
    '''
    if eta_roof < eta_roof_thr:
        return eta_insscr*( U_thscr*f2_vent_side \
            +(1 - U_thscr)*f_vent_roof_side*eta_side ) + 0.5*f_leakage
    else:
        return eta_insscr*f2_vent_side + 0.5*f_leakage

def f2_vent_side(C_d, U_side, A_side, v_wind, A_flr, C_w):
    '''
    Calculate the ventilation rate for sidewalls ventilation only

    :param C_d: Discharge coefficent
    :param U_side: Control of the side ventilators
    :param A_side: Maximum sidewall ventilation area
    :param v_wind: The wind speed
    :param C_w: Global wind pressure coefficent
    :param A_flr: Greenhouse floor area
    '''
    return ((C_d*U_side*A_side*v_wind)/(2*A_flr))*np.sqrt(C_w)

def f_vent_forced(eta_insscr, U_vent_forced, phi_vent_forced, A_flr):
    '''
    Calculate the forced ventilation

    :param eta_insscr: The ventilation rate's reduce factor created by insect screens
    :param U_vent_forced: The control of the forced ventilation
    :param phi_vent_forced: The air flow capacity of the forced ventilation system
    :param A_flr: Greenhouse floor area
    '''
    return (eta_insscr*U_vent_forced*phi_vent_forced)/A_flr
#------------------------------------------------------------------


def MC_top_out(f_vent_roof, CO2_top, CO2_out):
    '''
    Calculate the CO2 exchange between the top compartment air and the outside air

    :param f_vent_roof: The ventilation rate through roof openings
    :param CO2_top: CO2-concentration of the top compartment air
    :param CO2_out: CO2-concentration of the outside air
    '''
    return f_vent_roof*(CO2_top - CO2_out)

def f_vent_roof(eta_insscr, f2_vent_roof, f_leakage, U_thscr, f_vent_roof_side,
        eta_roof, eta_roof_thr):
    '''
    Calculate the ventilation rate through roof openings

    :param eta_insscr: The ventilation rate's reduce factor created by insect screens
    :param f2_vent_roof: The ventilation rate for roof vent ventilation only
    :param f_leakage: The leakage rate
    :param U_thscr: Control of the thermal screen
    :param f_vent_roof_side: The ventilation rate through both the roof and the side vents
    :param eta_roof: The ratio between the roof vents area and total ventilation area
    :param eta_roof_thr: The threshold value for which there is no chimney effect
    '''
    if eta_roof < eta_roof_thr:
        return eta_insscr*( U_thscr*f2_vent_roof \
            +(1 - U_thscr)*f_vent_roof_side*eta_roof ) + 0.5*f_leakage
    else:
        return eta_insscr*f2_vent_roof + 0.5*f_leakage

def f2_vent_roof(C_d, U_roof, A_roof, A_flr, g, h_roof, T_air, T_out, T_mean_air, C_w, v_wind):
    '''
    Calculate the ventilation rate due to roof ventilation

    :param C_d: Discharge coefficent
    :param U_roof: Control of the aperture of the roof vent
    :param A_roof: Maximum roof ventilation area
    :param A_flr: Greenhouse floor area
    :param g: Gravitational accleration
    :param h_roof: The vertical dimension of a single ventilation opening
    :param T_air: Greenhouse air temperature
    :param T_out: Outdoor temperature
    :param T_mean_air: Mean temperature of the greenhouse air and the outside air
    :param C_w: Global wind pressure coefficent
    :param v_wind: The wind speed
    '''
    return ( (C_d*U_roof*A_roof)/(2*A_flr) ) \
        *( (g*h_roof*(T_air - T_out))/(2*(T_mean_air + 273.15)) + C_w*v_wind**2 )**(1/2)

def eta_roof(U_roof, U_side, A_roof, A_side):
    '''
    Calculate the ratio between the roof vents area and total ventilation area
    :param U_roof: Control of the aperture of the roof vent
    :param A_roof: Maximum roof ventilation area
    :param A_side: Maximum sidewall ventilation area
    '''
    # if U_roof > 0.5:
    #     U_roof = 0.5
    # if U_roof < 0:
    #     U_roof = 0
    # return (A_roof*U_roof)/(A_roof + A_side)

    if (A_roof*U_roof + A_side*U_side == 0):
        return 0
    return (A_roof*U_roof)/(A_roof*U_roof + U_side*A_side)

def eta_side(U_roof, U_side, A_roof, A_side):
    '''
    Calculate the ratio between the roof vents area and total ventilation area
    :param U_roof: Control of the aperture of the roof vent
    :param A_roof: Maximum roof ventilation area
    :param A_side: Maximum sidewall ventilation area
    '''
    # if U_side > 0.5:
    #     U_side = 0.5
    # if U_side < 0:
    #     U_side = 0
    # return (A_side*U_side)/(A_roof + A_side)

    if (A_roof*U_roof + A_side*U_side == 0):
        return 0
    return (A_side*U_side)/(A_roof*U_roof + U_side*A_side)
#------------------------------------------------------------------


def MC_air_can(M_CH2O, h_Cbuf, P, R_P):
    '''
    Calculate the CO2 flux from the air to the canopy

    :param M_CH2O: The molar mass of CH2O
    :param h_Cbuf: The inhibition of the photosynthesis rate by saturation of the leaves
    with carbonhydrates
    :param P: The photonsynthesis rate
    :param R: The photorespiration during the photosynthesis process
    '''
    return M_CH2O*h_Cbuf*(P - R_P)

def h_Cbuf(C_buf, C_max_buf):
    '''
    Calculate the inhibition of the photosynthesis rate by saturation of the leaves
    with carbonhydrates

    :param C_buf: The buffer capacity
    :param C_max_buf: The maximum buffer capacity
    '''
    if C_buf > C_max_buf:
        return 0
    else:
        return 1

def P(J, CO2_stom, Gamma):
    '''
    Calculate the canopy photosynthesis rate

    :param J: The electron transport rate
    :param CO2_stom: CO2-concentration in the stomata
    :param Gamma: CO2 compensation point
    '''
    return (J*(CO2_stom - Gamma))/(4*(CO2_stom + 2*Gamma))

def R_P(P, Gamma, CO2_stom):
    '''
    Calculate the hotorespiration during photosynthesis processes

    :param P: the canopy photosynthesis rate
    :param Gamma: CO2 compensation point
    :param CO2_stom: CO2-concentration in the stomata
    '''
    return P*(Gamma/CO2_stom)

def J(J_POT, alpha, PAR_can, Theta):
    '''
    Calculate the electron transport rate
    :param J_POT: The potential electron transport rate
    :param alpha: The conversions factor from photons to electrons
    :param PAR_can: The absorbed PAR
    :param Theta: The degree of curvature of electron transport rate
    '''
    return (J_POT + alpha*PAR_can - \
        np.sqrt( (J_POT + alpha*PAR_can)**2 - 4*Theta*J_POT \
            *alpha*PAR_can ))/ (2*Theta)

def J_POT(J_MAX_25can, Ej, T_can, T_25, R, S, H):
    '''
    Calculate the potential electron transport rate

    :param J_MAX_25can: The maximum rate of electron transport at 25oC for the canopy
    :param Ej: The activation energy for J_POT
    :param T_can: Canopy temperature (K)
    :param T_25: The reference temperature at 25oC (K)
    :param R: The molar gas constant
    :param S: The entropy term
    :param H: The deactivation energy
    '''
    return J_MAX_25can*np.exp(Ej*( (T_can + 273.15) -T_25)/(R*(10**-3)* (T_can + 273.15) *T_25)) \
        * (1 + np.exp((S*T_25 - H)/(R*(10**-3)*T_25))) \
        /(1 + np.exp((S* (T_can + 273.15) - H)/(R*(10**-3)* (T_can + 273.15) )))

def J_MAX_25can(LAI, J_MAX_25leaf):
    '''
    Calculate the maximum rate of electron transport at 25oC for the canopy

    :param LAI: leaf area index
    :param J_MAX_25leaf: The maximum rate of electron transport for the
    leaf at 25oC
    '''
    return LAI*J_MAX_25leaf

def CO2_stom(eta_CO2_air_stom, CO2_air):
    '''
    Calculate the CO2-concentration in the stomata

    :param eta_CO2_air_stom: conversion factor
    :param CO2_air: CO2-concentration of the greenhouse air
    '''
    return eta_CO2_air_stom*CO2_air

def Gamma(J_MAX_25can, J_MAX_25leaf, c_Gamma, T_can):
    '''
    Calculate the CO2 compensation point

    :param J_MAX_25can: The maximum rate of electron transport at 25oC for the canopy
    :param J_MAX_25leaf: The maximum rate of electron transport for the
    :param c_Gamma: the effect of canopy temperature on the CO2 compensation point
    :param T_can: The canopy temperature
    '''
    return (J_MAX_25leaf/J_MAX_25can)*c_Gamma*(T_can + 273.15) \
        + 20*c_Gamma*(1 - (J_MAX_25leaf/J_MAX_25can))
#------------------------------------------------------------------

