'''
===========================================
Python script for Math modelling assignment
===========================================
'''

''' Dependencies '''

import math
import numpy as np
'''
TODO: Find the value for the parameters, and write the dx function
'''

''' ------------------------ '''
''' Time-variant params      '''
''' ------------------------ '''
T_air = 20
T_out = 10
T_mean_air = 15
T_top = 10
T_mech_cool = 0
T_thscr = 20
T_cov_in = 9
T_can = 10
valve =0.15
U_fog = valve
U_mech_cool = valve
U_pad = valve
U_blow = valve
U_side = valve
U_vent_forced = valve
U_ThScr = 0.97
U_roof = 1
U_shscr = valve
U_thscr = U_ThScr
''' ------------------------ '''
''' Time-invariant params    '''
''' ------------------------ '''
v_wind = 0
eta_heatCO2 = 0.057
A_flr = 14000
A_roof = 0.1*A_flr
A_side = 0
h_vent = 0.004
P_blow = 1
g = 9.81
phi_extCO2 = 7.2*10**4
K_thscr = 0.05*10**-3
M_air = 28.96
h_elevator = 0
rho_air0 = 1.2
C_d_Gh = 0.75
C_w_Gh = 0.09
zeta_insscr = 1
c_leakage = 10**-4
eta_side = 0
M_CH2O = 30*10**-3
C_bufmax = 20*10**-3
H_d = 22*10**4
H_a = 37*10**3
R = 8.314*1000
S = 710
m = 0.1
tau_Gh = 0.78
eta_glob_par = 2.3
P_MLT = 50
L_05 = 406
eta_heatVap = 4.43*10**-8
phi_fog = 0
eta_pad = 0
phi_pad = 0
M_water = 18
x_pad = 0
x_out = 0
s_mv_12 = -0.1
delta_H = 2.45*10**6
c_HEC_in = 1.86
c_p_air = 1000
gamma = 65.8
r_b = 275
r_s_min = 200
I_glob = 179.4
R_can = 0
c_evap1 = 4.3
c_evap2 = 0.54
c_night_evap3 = 1.1*10**(-11)
c_night_evap4 = 5.2*10**(-6)
c_day_evap3 = 6.1*10**(-7)
c_day_evap4 = 5.2*10**(-6)
eta_mmg_ppm = 1
srs = -1
R_can_SP = 5
eta_roof_thr = 0.9
zeta_insscr = 0.9
rho_air0 = 1.2
h_elevation = 0
LAI = 2.5
phi_pad = 0
eta_shscrCd = 0
eta_shscrCw = 0
phi_vent_forced = 0
A_cov = 18000
COP_mech_cool = 0
P_mech_cool = 0
h_air = 3.8
h_top = 0.4
h_thscr = 0.35*10**(-3)
h_side_roof = 0
h_elevation_top = 3.8
h_elevation_mean_air = 4.2
VP_air = 1
VP_top = 2
VP_out = 2100
CO2_air = 417


def rho(h):
    '''
    Compute the saturated Vapour Pressure

    :param h: Altitude of place need to be computed 
    '''
    return 101325 *(1 - 2.25577* 10**(-5)* h)** 5.25588
def rho_air():
    '''
    Calculte the density of the greenhouse air

    :param rho_air0: The density of the air at sea level
    :param g: Gravitational accleration
    :param M_air: The molar mass of air
    :param h_elevation: The altitude of the greenhouse above sea level
    :param T_air: Greenhouse air temperature
    :param R: The molar gas constant
    '''
    return M_air*rho(h_elevation)/((273.15+T_air)*R)
    # return rho_air0*np.exp( (g*M_air*h_elevation)/((273.15+T_air)*R) )

def rho_out():
    '''
    Calculate the density of the outdoor air

    :param rho_air0: The density of the air at sea level
    :param g: Gravitational accleration
    :param M_air: The molar mass of air
    :param h_elevation: The altitude of the greenhouse above sea level
    :param T_out: Outdoor temperature
    :param R: The molar gas constant
    '''
    return M_air*rho(h_elevation)/((273.15+T_out)*R)
    # return rho_air0*np.exp( (g*M_air*h_elevation)/((273.15+T_out)*R) )

def rho_top():
    '''
    Calculate the density of the outdoor air

    :param rho_air0: The density of the air at sea level
    :param g: Gravitational accleration
    :param M_air: The molar mass of air
    :param h_elevation_top: The altitude of the top compartment above sea level
    :param T_top: Top compartment temperature
    :param R: The molar gas constant
    '''
    return M_air*rho(h_elevation)/((273.15+T_top)*R)
    # return rho_air0*np.exp((g*M_air*h_elevation_top)/((273.15+T_top)*R))

def rho_mean_air():
    '''
    Calculate the mean density of the top and below compartment of the greenhouse air

    :param rho_air0: The density of the air at sea level
    :param g: Gravitational accleration
    :param M_air: The molar mass of air
    :param h_elevation_mean_air: The mean altitude of both compartment 
    :param T_mean_air: Mean temperature 
    :param R: The molar gas constant
    '''
    return M_air*rho(h_elevation)/((273.15+T_mean_air)*R)
    # return rho_air0*np.exp((g*M_air*h_elevation_mean_air)/((273.15+T_mean_air)*R))

''' (7) '''
def f_thscr():
    '''
    Calculate the air flux rate through the thermal screen 

    :param U_thscr: Control of the thermal screen
    :param K_thscr: The thermal screen flux coefficent
    :param T_air: Greenhouse air temperature
    :param T_top: Greenhouse top compartment air temperature
    :param g: Gravitational accleration
    :param rho_mean_air: The mean density of the greenhouse and the top compartment air
    :param rho_air: The density of greenhouse air
    :param rho_top: The density of top compartment air
    '''
    return U_thscr*K_thscr*(abs(T_air - T_top))**(2/3) + \
        (1 - U_thscr)/rho_mean_air()*( (g*(1-U_thscr))/(2*rho_mean_air())*abs(rho_air() - rho_top()) )**(1/2)


def C_d():
    '''
    Calculate the discharge coefficent

    :param C_d_Gh: the discharge coefficient determined for a greenhouse 
    without an external shading screen
    :param eta_shscrCd: parameter that determines the effect of the 
    shading screen on the discharge coefficient
    :param U_shscr: Control of the external shading screen
    '''
    return C_d_Gh*(1 - eta_shscrCd*U_shscr)

def C_w():
    '''
    Calculate the global wind pressure coefficent

    :param C_w_Gh: the global wind pressure coefficient for a greenhouse
    without an external shading screen
    :param eta_shscrCw: parameter that determines the effect of the 
    shading screen on the global wind pressure coefficient
    :param U_shscr: Control of the external shading screen
    '''
    return C_w_Gh*(1 - eta_shscrCw*U_shscr)

''' (10) '''
def f_vent_roof_side():
    '''
    Calculate the ventilation rate through both the roof and the side vents

    :param C_d: Discharge coefficent
    :param A_flr: Greenhouse floor area
    :param U_roof: Control of the aperture of the roof vent
    :param U_side: Control of the side ventilators
    :param A_roof: Maximum roof ventilation area
    :param A_side: Maximum sidewall ventilation area
    :param g: Gravitational accleration
    :param h_side_roof: Vertical distance between mid-points of side wall and 
    roof ventilation openings
    :param T_air: Greenhouse air temperature
    :param T_out: Outdoor temperature
    :param T_mean_air: Mean temperature of the greenhouse air and the outside air 
    :param C_w: Global wind pressure coefficent
    :param v_wind: The wind speed
    '''
    if (U_roof*A_roof==0):
        return 0
    return (C_d()/A_flr)* \
        ( ((U_roof**2)*(U_side**2)*(A_roof**2)*(A_side**2))/ \
        ((U_roof**2)*(A_roof**2)+(U_side**2)*(A_side**2)) \
        *(2*g*h_side_roof*(T_air - T_out))/(T_mean_air+273.15) \
        + (((U_roof*A_roof + U_side*A_side)/2)**2)*C_w()*v_wind**2  )**(1/2)
''' (11) '''
def eta_insscr():
    '''
    Calculate the ventilation rate's reduce factor created by insect screens

    :param zeta_insscr: The screen porosity
    '''
    return zeta_insscr*(2 - zeta_insscr)
''' (12) '''
def f_leakage():
    '''
    Calculate the leakage rate

    :param c_leakage: The leakage coefficent
    :param v_wind: The wind speed
    '''
    if v_wind < 0.25:
        return 0.25*c_leakage
    else:
        return v_wind*c_leakage
''' (13) '''
def f_vent_side():
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
    if eta_roof() < eta_roof_thr:
        return eta_insscr()*( U_thscr*f2_vent_side() \
            +(1 - U_thscr)*f_vent_roof_side()*eta_side()) + 0.5*f_leakage()
    else:
        return eta_insscr()*f2_vent_side() + 0.5*f_leakage()
def f2_vent_side():
    '''
    Calculate the ventilation rate for sidewalls ventilation only

    :param C_d: Discharge coefficent
    :param U_side: Control of the side ventilators
    :param A_side: Maximum sidewall ventilation area
    :param v_wind: The wind speed
    :param C_w: Global wind pressure coefficent
    :param A_flr: Greenhouse floor area
    '''
    return ((C_d()*U_side*A_side*v_wind)/(2*A_flr))*np.sqrt(C_w())

''' (14) '''
def f_vent_forced():
    '''
    Calculate the forced ventilation

    :param eta_insscr: The ventilation rate's reduce factor created by insect screens
    :param U_vent_forced: The control of the forced ventilation
    :param phi_vent_forced: The air flow capacity of the forced ventilation system
    :param A_flr: Greenhouse floor area
    '''
    return (eta_insscr()*U_vent_forced*phi_vent_forced)/A_flr

def eta_roof():
    '''
    Calculate the ratio between the roof vents area and total ventilation area
    :param U_roof: Control of the aperture of the roof vent
    :param A_roof: Maximum roof ventilation area
    :param A_side: Maximum sidewall ventilation area
    '''
    if (A_roof + A_side == 0):
        return 0
    return (A_roof*U_roof)/(A_roof + A_side) 

def eta_side():
    '''
    Calculate the ratio between the roof vents area and total ventilation area
    :param U_roof: Control of the aperture of the roof vent
    :param A_roof: Maximum roof ventilation area
    :param A_side: Maximum sidewall ventilation area
    '''
    if (A_roof + A_side == 0):
        return 0
    return (A_side*U_side)/(A_roof+ A_side)

''' (16) '''
def f_vent_roof():
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
    if eta_roof() < eta_roof_thr:
        return eta_insscr()*( U_thscr*f2_vent_roof() \
            +(1 - U_thscr)*f_vent_roof_side()*eta_roof() ) + 0.5*f_leakage()
    else:
        return eta_insscr()*f2_vent_roof() + 0.5*f_leakage()
''' (17) '''
def f2_vent_roof():
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
    return ( (C_d()*U_roof*A_roof)/(2*A_flr) ) \
        *( (g*h_vent*(T_air - T_out))/(2*(T_mean_air+273.15)) + C_w()*v_wind**2 )**(1/2)



''' (6) '''
def MV_blow_air():
    '''
    Calculate the flow of vapour from blower to the greenhouse air

    :param eta_heatvap: Amount of vapour which is released when 1 Joule of sensible energy is produced by the heat blower 
    :param U_blow: 		The control value of the direct air heater
    :param P_blow: 		The heat capacity of the direct air heater
    :param A_flr: 		Greenhouse floor area
    '''
    return (eta_heatVap*U_blow*P_blow)/A_flr

''' (7) '''
def MV_fog_air():
	'''
	Caculate amount of vapour supplied from fogging system
    :param U_fog: 	The control value of the direct fogging system
    :param phi_fog: fogging system maximum supply
    :param A_flr: 	Greenhouse floor area
	'''
	return (U_fog*phi_fog)/A_flr

''' (8) '''
def MV_pad_air():
	'''
	Caculate the ability of air crossing the pad and fan system
	:param rho_air:	The density of air
	:param U_pad:	The control value of pad and fan system
	:param phi_pad: Pad and fan system maximum supply
	:param A_flr: 	Greenhouse floor area
	:param eta_pad: Efficiency of the pad and fan system
	:param x_pad:	Water vapour content in pad and fan system
	:param x_out:	Outdoor water vapour content 
	'''
	return (rho_air()*U_pad*phi_pad*(eta_pad*(x_pad - x_out)+x_out))/A_flr

''' (9) '''
def MV_air_out_pad():
	'''
	Caculate amount of vapour flux loss out
	:param U_pad:	The control value of pad and fan system
	:param phi_pad:	Pad and fan system maximum supply
	:param M_water: Molar mass of water
	:param VP_air:	Air vapour pressure 
	:param T_air:	Air temperature
	:param A_flr:	Greenhouse floor area
	:param R:		FIR flux density
	'''
	return (U_pad*phi_pad*M_water*VP_air)/(A_flr*R*(T_air+273.15))

''' (10)'''
def MV_air_mech():
	'''
	Caculate The amount of vapour in the air is also affected by system mechanical cool
	:param U_mech_cool:	  Control of the mechanical cooling
	:param COP_mech_cool: Coefficient of performance of the mechanical cooling system
	:param P_mech_cool:   Electrical capacity of the mechanical cooling system
	:param A_flr:		  Greenhouse floor area
	:param T_air:		  Air temperature
	:param T_mech_cool:	  The temperature of the cool surface of the mechanical cooling system
	:param VP_air:		  Air vapour pressure
	:param VP_mech_cool:  The mechanical cooling system vapour pressure
	:param xic_ma:		  Stefan Boltzmann constant
	:param H:		      Sensible heat flux density
	'''
	return (6.04*pow(10,-9)*((U_mech_cool*COP_mech_cool*P_mech_cool)/A_flr)*(VP_air-VP_mech()))/(T_air-T_mech_cool+6.4*pow(10,-9)*delta_H*(VP_air-VP_mech()))
'''(11)'''
def MV_air_thscr():
    '''
    Caculate the amount of vapour lost due to condensation forming a heat shield
    :param HEC_air_thscr: Heat exchange coefficient 
    :param VP_air:		  Air vapour pressure 
    :param VP_mech:		  Mechanical cooling vapour pressure
    '''
   
    return 6.04*pow(10,-9)*HEC_air_thscr()*(VP_air-VP_ThScr())/(1+np.exp(s_mv_12*(VP_air-VP_ThScr())))

'''(12)'''
def HEC_air_thscr():
	'''
	Caculate heat exchange coefficient
	:param U_thscr:	Control of the thermal screen 
	:param T_air:	Air temperature
	:param T_thscr:	Thermal screen temperature
	'''
	return 1.7*U_thscr*abs(T_air-T_thscr)

'''(13)'''
def MV_top_cov_in():
    '''
    Caculate mass vapour flux cover internal side
    :param VP_air:		  Air vapour pressure 
    :param VP_mech:		  Mechanical cooling vapour pressure
    '''
    return 6.04*pow(10,-9)*HEC_top_cov_in()*(VP_air-VP_Top_cov_in())/(1+np.exp(s_mv_12*(VP_air-VP_Top_cov_in())))

'''(14)'''
def HEC_top_cov_in():
	'''
	Caculate heat exchange coefficient top cover internal side
	:param c_HEC_in: Coefficient of heat exchange between the cover and the air in the environment outside
	:param T_top:	 Tempature of compartment above the thermal screen
	:param T_cov_in: Tempature of cover internal side
	:param A_cov: 	 Greenhouse cover area
	:param A_flr: 	 Greenhouse floor area
	'''
	return (c_HEC_in*pow(T_top - T_cov_in,0.33)*A_cov)/A_flr

'''(15)'''
def MV_air_top():
	'''
	Caculate the amount of vapour going from air to top
	:param M_water:	Molar mass of water
	:param R:		FIR flux density
	:param f_thscr: Caculate form f_thscr funtion
	:param VP_air: 	Air vapour pressure 
	:param T_air: 	Air temperature
	:param VP_top: 	Compartment above the thermal screen vapour pressure
	:param T_top: 	Tempature of compartment above the thermal screen
	'''
	return (M_water*f_thscr()*((VP_air)/(T_air+273.15) - (VP_top)/(T_top+273.15)))/R

'''(16)'''
def MV_air_out():
	'''
	Caculate the amount of vapour going from air to outside
	:param M_water:		  Molar mass of water
	:param R:			  FIR flux density
	:param f_vent_side:   Caculate form f_vent_side funtion
	:param f_vent_forced: Caculate form f_vent_forced funtion
	:param VP_air: 		  Air vapour pressure 
	:param T_air: 		  Air temperature
	:param VP_out: 		  Outside vapour pressure
	:param T_out: 		  Tempature of outside
	'''	
	return (M_water*(f_vent_side()+f_vent_forced())*((VP_air)/(T_air+273.15) - (VP_out)/(T_out+273.15)))/R

'''(17)'''
def MV_top_out():
	'''
	Caculate the amount of vapour going from top to outside
	:param M_water:		  Molar mass of water
	:param R:			  FIR flux density
	:param f_vent_roof:   Caculate form f_vent_roof funtion
	:param VP_top: 		  Top vapour pressure 
	:param T_top: 		  Top temperature
	:param VP_out: 		  Outside vapour pressure
	:param T_out: 		  Tempature of outside
	'''
	return (M_water*f_vent_roof()*((VP_top)/(T_top+273.15) - (VP_out)/(T_out+273.15)))/R


'''(18)'''
def VEC_can_air():
	'''
	Caculate vapour exchange coefficient
	:param rho_air: The density of greenhouse air
	:param c_p_air: Heat capacity of air
	:param LAI:		The leaf area index
	:param delta_h: Latent heat of evaporation
	:param gamma:		Psychrometric constant
	:param r_b:		Boundary layer resistance of the canopy for vapour transport
	:param r_s:		The canopy resistance for transpiration
	'''
	return (2*rho_air()*c_p_air*LAI)/(delta_H*gamma*(r_b+r_s()))
def MV_can_air():
	'''
	Caculate the evaporation of the foliage
	:param VEC_can_air: Caculate from VEC_can_air funtion
	:param VP_can: 		Canopy vapour pressure
	:param VP_air:		Air vapour pressure 
	'''
	return VEC_can_air()*(VP_can()-VP_air)

def r_s():
    '''
	Caculate the transpiration resistance
	:param r_s_min:    Minimum canopy resistance 
	'''
    # return 82
    return r_s_min*rf_R_can()*rf_CO2()*rf_VP()
def rf_R_can():
    '''
	Caculate the canopy resistance factor
	:param R_can:       Radiation above the canopy
	:param c_evap1: 	Observed param 1
	:param c_evap2:		Observed param 1 
	'''
    return ((R_can+c_evap1)/(R_can+c_evap2))
def rf_CO2():
    '''
	Caculate the resistance factor of Co2 in the lower compartment
	:param c_evap3: 	Observed param 3
	:param CO2_air:		CO2 amount in the lower compartment 
	'''
    temp = 1+c_evap3()*(CO2_air-200)**2
    return temp
    if temp <= 1.5:
        return temp
    return 1.5
    # return temp
def rf_VP():
    '''
	Caculate the resistance factor of vapour pressire in the lower compartment
	:param c_evap3: 	Observed param 4
	:param VP_air:		Vapour pressure in the lower compartment 
    :param VP_can:		Vapour pressure at the canopy 
	'''
    temp = 1+c_evap4()*(VP_can()-VP_air)**2
    return temp
    if temp <=5.8:
        return temp
    return 5.8
    # return temp
def S_rs():
    '''
	Caculate the switch function
	:param s_r_s: 	    Slope of the function
	:param R_can:		Radiation above canopy 
    :param R_can_SP:	Radiation abobe canopy at night 
	'''
    return 1/(1+np.exp(srs*(R_can-R_can_SP)))

def c_evap3():
    '''
	Caculate the resistance factor of vapour pressire in the lower compartment
	:param c_night_evap3: 	c_evap at night
	:param c_day_evap3:		c_evap at day
	'''
    return c_night_evap3*(1-S_rs())+c_day_evap3*S_rs()
def c_evap4():
    '''
	Caculate the resistance factor of vapour pressire in the lower compartment
	:param c_night_evap4: 	c_evap at night
	:param c_day_evap4:		c_evap at day
	'''
    return c_night_evap4*(1-S_rs())+c_day_evap4*S_rs()

def cap_VP_air():
    return (M_water*h_air)/(R*(T_air+273.15))
def cap_VP_top():
    return (M_water*h_top)/(R*(T_top+273.15))
def VP_sat(t):
    '''
	Caculate saturated vapour pressure
	:param t: Temperature		
	'''
    return 610.78 * np.exp( t / ( t + 238.3 ) * 17.2694 )

def VP_ThScr():
    '''
	Caculate saturated vapour pressire on the thermal screen
	:param T_thscr: Temperature of thermal screen		
	'''
    return VP_sat(T_thscr)
    # return (-274.36+877.52*math.exp(0.0545*T_thscr))
def VP_mech():
    '''
	Caculate saturated vapour pressire on the mechanical cooling surface
	:param T_mech_cool: Temperature of mechanical cooling surface	
	'''
    return VP_sat(T_mech_cool)
    # return (-274.36+877.52*math.exp(0.0545*T_mech_cool))
def VP_Top_cov_in():
    '''
	Caculate saturated vapour pressire on the top covering layer
	:param T_cov_in: Temperature of thermal screen		
	'''
    return VP_sat(T_cov_in)
    # return (-274.36+877.52*math.exp(0.0545*T_cov_in))
def VP_can():
    '''
	Caculate saturated vapour pressire in the canopy
	:param T_can: Temperature of canopy		
	'''
    return VP_sat(T_can)
    # return (-274.36+877.52*math.exp(0.0545*T_can))
def VP_Out():
    '''
	Caculate saturated vapour pressire outside
	:param T_out: Temperature of canopy		
	'''
    return VP_sat(T_out)
    # return (-274.36+877.52*math.exp(0.0545*T_out))
def VP_Top():
    '''
	Caculate saturated vapour pressire on the top compartment
	:param T_top: Temperature of top compartment
	'''
    return VP_sat(T_top) 
    # return (-274.36+877.52*math.exp(0.0545*T_top))
def VP_Air():
    '''
	Caculate saturated vapour pressire in the lower compartment
	:param T_air: Temperature of lower compartment		
	'''
    return VP_sat(T_air)
    # return (-274.36+877.52*math.exp(0.0545*T_air))
