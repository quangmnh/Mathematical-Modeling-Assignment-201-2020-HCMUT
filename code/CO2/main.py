import init
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.stdout = open("log.txt", "w")

list_CO2_air_euler = []
list_CO2_top_euler = []
list_CO2_air_rk4 = []
list_CO2_top_rk4 = []
list_real_CO2_air_data = []
list_err_euler = []
list_err_rk4 = []

list_MC_ext_air = []
list_MC_air_top_Air = []
list_MC_air_out = []
list_MC_air_can = []

list_MC_air_top_Top = []
list_MC_top_out = []

def randomU():
    init.U_blow = np.random.uniform(0, 1)
    init.U_pad = np.random.uniform(0, 1)
    init.U_shscr = np.random.uniform(0, 1)
    init.PAR_can = np.random.normal(100, 50)

def f_func(CO2_air, CO2_top, flag):

    MC_blow_air = init.MC_blow_air(
        init.eta_heatCO2,
        init.U_blow,
        init.P_blow,
        init.A_flr
    )

    MC_ext_air = init.MC_ext_air(
        init.U_extCO2,
        init.phi_extCO2,
        init.A_flr
    )

    MC_pad_air = init.MC_pad_air(
        init.U_pad,
        init.phi_pad,
        init.A_flr,
        init.CO2_out,
        CO2_air
    )

    MC_air_can = init.MC_air_can(
        init.M_CH2O,
        init.h_Cbuf(init.C_buf, init.C_max_buf),
        init.P(
            init.J(
                init.J_POT(
                    init.J_MAX_25can(init.LAI, init.J_MAX_25leaf),
                    init.Ej,
                    init.T_can,
                    init.T_25,
                    init.R,
                    init.S,
                    init.H
                ),
                init.alpha,
                init.PAR_can,
                init.Theta
            ),
            init.CO2_stom(
                init.eta_CO2_air_stom,
                CO2_air
            ),
            init.Gamma(
                init.J_MAX_25can(init.LAI, init.J_MAX_25leaf),
                init.J_MAX_25leaf,
                init.c_Gamma,
                init.T_can
            )
        ),
        init.R_P(
            init.P(
                init.J(
                    init.J_POT(
                        init.J_MAX_25can(
                            init.LAI,
                            init.J_MAX_25leaf
                        ),
                        init.Ej,
                        init.T_can,
                        init.T_25,
                        init.R,
                        init.S,
                        init.H
                    ),
                    init.alpha,
                    init.PAR_can,
                    init.Theta
                ),
                init.CO2_stom(
                    init.eta_CO2_air_stom,
                    CO2_air
                ),
                init.Gamma(
                    init.J_MAX_25can(init.LAI, init.J_MAX_25leaf),
                    init.J_MAX_25leaf,
                    init.c_Gamma,
                    init.T_can
                )
            ),
            init.Gamma(
                init.J_MAX_25can(init.LAI, init.J_MAX_25leaf),
                init.J_MAX_25leaf,
                init.c_Gamma,
                init.T_can
            ),
            init.CO2_stom(
                init.eta_CO2_air_stom,
                CO2_air
            )
        )
    )

    MC_air_top = init.MC_air_top(
        init.f_thscr(
            init.U_thscr,
            init.K_thscr,
            init.T_air,
            init.T_top,
            init.g,
            init.rho_mean_air(
                init.rho_air0,
                init.g,
                init.M_air,
                init.h_elevation_mean_air,
                init.T_air,
                init.T_top,
                init.R
            ),
            init.rho_air(
                init.rho_air0,
                init.g,
                init.M_air,
                init.h_elevation,
                init.T_air,
                init.R
            ),
            init.rho_top(
                init.rho_air0,
                init.g,
                init.M_air,
                init.h_elevation_top,
                init.T_top,
                init.R
            )
        ),
        CO2_air,
        CO2_top
    )

    MC_air_out = init.MC_air_out(
        init.f_vent_side(
            init.eta_insscr(init.zeta_insscr),
            init.f2_vent_side(
                init.C_d(
                    init.C_d_Gh,
                    init.eta_shscrCd,
                    init.U_shscr
                ),
                init.U_side,
                init.A_side,
                init.v_wind,
                init.A_flr,
                init.C_w(
                    init.C_w_Gh,
                    init.eta_shscrCw,
                    init.U_shscr
                )
            ),
            init.f_leakage(
                init.c_leakage,
                init.v_wind
            ),
            init.U_thscr,
            init.f_vent_roof_side(
                init.C_d(
                    init.C_d_Gh,
                    init.eta_shscrCd,
                    init.U_shscr
                ),
                init.A_flr,
                init.U_roof,
                init.U_side,
                init.A_roof,
                init.A_side,
                init.g,
                init.h_side_roof,
                init.T_air,
                init.T_out,
                init.T_mean_air(init.T_air, init.T_out),
                init.C_w(
                    init.C_w_Gh,
                    init.eta_shscrCw,
                    init.U_shscr
                ),
                init.v_wind
            ),
            init.eta_side(
                init.U_roof,
                init.U_side,
                init.A_roof,
                init.A_side
            ),
            init.eta_roof(
                init.U_roof,
                init.U_side,
                init.A_roof,
                init.A_side
            ),
            init.eta_roof_thr
        ),
        init.f_vent_forced(
            init.eta_insscr(init.zeta_insscr),
            init.U_vent_forced,
            init.phi_vent_forced,
            init.A_flr
        ),
        CO2_air,
        init.CO2_out
    )

    if flag:
        # print("MC_ext_air= ", MC_ext_air)
        # print("MC_air_top= ", MC_air_top)
        # print("MC_air_out= ", MC_air_out)
        # print("MC_air_can= ", MC_air_can)
        # print("Result= ", (MC_blow_air + MC_ext_air + MC_pad_air - MC_air_can - MC_air_top - MC_air_out) / init.cap_CO2_air)
        # print("--------------------------")

        list_MC_air_can.append(MC_air_can)
        list_MC_air_top_Air.append(MC_air_top)
        list_MC_air_out.append(MC_air_out)
        list_MC_ext_air.append(MC_ext_air)

    return (MC_blow_air + MC_ext_air + MC_pad_air - MC_air_can - MC_air_top - MC_air_out) / init.cap_CO2_air


def g_func(CO2_air, CO2_top, flag):
    MC_air_top = init.MC_air_top(
        init.f_thscr(
            init.U_thscr,
            init.K_thscr,
            init.T_air,
            init.T_top,
            init.g,
            init.rho_mean_air(
                init.rho_air0,
                init.g,
                init.M_air,
                init.h_elevation_mean_air,
                init.T_air,
                init.T_top,
                init.R
            ),
            init.rho_air(
                init.rho_air0,
                init.g,
                init.M_air,
                init.h_elevation,
                init.T_air,
                init.R
            ),
            init.rho_top(
                init.rho_air0,
                init.g,
                init.M_air,
                init.h_elevation_top,
                init.T_top,
                init.R
            )
        ),
        CO2_air,
        CO2_top
    )

    MC_top_out = init.MC_top_out(
        init.f_vent_roof(
            init.eta_insscr(init.zeta_insscr),
            init.f2_vent_roof(
                init.C_d(
                    init.C_d_Gh,
                    init.eta_shscrCd,
                    init.U_shscr
                ),
                init.U_roof,
                init.A_roof,
                init.A_flr,
                init.g,
                init.h_roof,
                init.T_air,
                init.T_out,
                init.T_mean_air(
                    init.T_air,
                    init.T_out
                ),
                init.C_w(
                    init.C_w_Gh,
                    init.eta_shscrCw,
                    init.U_shscr
                ),
                init.v_wind
            ),
            init.f_leakage(
                init.c_leakage,
                init.v_wind
            ),
            init.U_thscr,
            init.f_vent_roof_side(
                init.C_d(
                    init.C_d_Gh,
                    init.eta_shscrCd,
                    init.U_shscr
                ),
                init.A_flr,
                init.U_roof,
                init.U_side,
                init.A_roof,
                init.A_side,
                init.g,
                init.h_side_roof,
                init.T_air,
                init.T_out,
                init.T_mean_air(init.T_air, init.T_out),
                init.C_w(
                    init.C_w_Gh,
                    init.eta_shscrCw,
                    init.U_shscr
                ),
                init.v_wind
            ),
            init.eta_roof(
                init.U_roof,
                init.U_side,
                init.A_roof,
                init.A_side
            ),
            init.eta_roof_thr
        ),
        CO2_top,
        init.CO2_out
    )

    if flag:
        # print("MC_air_top= ", MC_air_top)
        # print("MC_top_out= ", MC_top_out)
        # print("-------------------------")
        list_MC_air_top_Top.append(MC_air_top)
        list_MC_top_out.append(MC_top_out)

    return (MC_air_top - MC_top_out) / init.cap_CO2_top


def dx(x0, y0, flag):
    '''
    Calculate the rate of change in CO2 concentration in greenhouse and top of greenhouse

    :param x0: CO2-concentration of the greenhouse air at time start t
    :param y0: CO2-concentration of the top greenhouse air at time start t
    :param flag: To be used for printing result to log file and saving to data for plotting
    '''
    CO2_air = f_func(x0, y0, flag)
    CO2_top = g_func(x0, y0, flag)
    return CO2_air, CO2_top


def euler(dx, x0, y0, h, flag):
    '''
    Calculate the approximate value of CO2-concentration of the greenhouse air
    and CO2-concentration of the top greenhouse air

    :param dx: the callable function of the ODE system
    :param x0: CO2-concentration of the greenhouse air at time start t
    :param y0: CO2-concentration of the top greenhouse air at time start t
    :param h: step approximation
    :param flag: To be used for printing result to log file and saving to data for plotting
    '''

    Kx, Ky = dx(x0, y0, flag)

    Kx = h * Kx
    Ky = h * Ky

    x0 = x0 + Kx
    y0 = y0 + Ky

    return x0, y0

def euler_loop(x0, y0, h, n):
    '''
        Calculate the approximate value of CO2-concentration of the greenhouse air
        and CO2-concentration of the top greenhouse air at t + h, t + 2h, ..., t + n*h

        :param x0: CO2-concentration of the greenhouse air at time start t
        :param y0: CO2-concentration of the top greenhouse air at time start t
        :param h: step approximation
        :param n: the number of times you want to get data (means you can still get data after
                  n*h (minutes))
    '''

    print("------EULER RESULT------")

    file_meteo = open("meteo.csv")
    file_meteo.readline()
    file_ghc = open("Greenhouse_climate.csv")
    file_ghc.readline()
    file_vip = open("vip.csv")
    file_vip.readline()

    CO2_air = x0
    CO2_top = y0
    step_next = 1
    step = 0.1
    n *= 300
    i = 0
    t0 = 0
    getData = False
    CO2_air_real_data = 427

    while i < n:
        if i // 300 >= step_next:
            step_next += 1

            line_meteo = file_meteo.readline()
            line_ghc = file_ghc.readline()
            line_vip = file_vip.readline()
            elements_vip = line_vip.split(',')

            while "NaN" in line_meteo or "NaN" in line_ghc or elements_vip[1] == "NaN":
                line_meteo = file_meteo.readline()
                line_ghc = file_ghc.readline()
                line_vip = file_vip.readline()
                elements_vip = line_vip.split(',')

            elements_vip = line_vip.split(',')
            elements_meteo = line_meteo.split(',')
            elements_ghc = line_ghc.split(',')

            init.v_wind = float(elements_meteo[10])
            init.T_out = float(elements_meteo[8])
            init.T_air = float(elements_ghc[9])
            init.U_thscr = float(elements_ghc[3]) / 100
            init.U_roof = (float(elements_ghc[10]) + float(elements_ghc[11])) / 200
            init.T_can = init.T_air + 15
            init.T_top = init.T_air + 1

            if elements_vip[0] == "NaN":
                init.U_extCO2 = np.random.uniform(0, 0.3)
            elif init.convert_mgm3_to_ppm(CO2_air) < float(elements_vip[0]):
                init.U_extCO2 = 1
            else:
                init.U_extCO2 = 0

            CO2_air_real_data = float(elements_ghc[2])

            getData = True

        if getData:
            randomU()
            getData = False

        CO2_air_step, CO2_top_step = euler(dx, CO2_air, CO2_top, step, getData)
        CO2_air = CO2_air_step
        CO2_top = CO2_top_step

        if t0 % (h * 60) == 0:

            print("CO2_real: ", CO2_air_real_data)
            list_real_CO2_air_data.append(CO2_air_real_data)

            err = abs(CO2_air_real_data - CO2_air)
            list_err_euler.append(init.convert_mgm3_to_ppm(err))

            list_CO2_air_euler.append(init.convert_mgm3_to_ppm(CO2_air))
            list_CO2_top_euler.append(init.convert_mgm3_to_ppm(CO2_top))
            print("CO2_air at t +", t0 // 60, "=", init.convert_mgm3_to_ppm(CO2_air))
            print("CO2_air at t +", t0 // 60, "=", init.convert_mgm3_to_ppm(CO2_top))
            print("------------------")

        i += 1
        t0 += 1

    file_meteo.close()
    file_ghc.close()
    file_vip.close()


def rk4(dx, x0, y0, h, flag):
    '''
        Calculate the approximate value of CO2-concentration of the greenhouse air
        and CO2-concentration of the top greenhouse air

        :param dx: the callable function of the ODE system
        :param x0: CO2-concentration of the greenhouse air at time start t
        :param y0: CO2-concentration of the top greenhouse air at time start t
        :param h: step approximation
        :param flag: To be used for printing result to log file and saving to data for plotting
        '''

    K1x, K1y = dx(x0, y0, flag)
    K1x = h * K1x
    K1y = h * K1y

    K2x, K2y = dx(x0 + K1x/2, y0 + K1y/2, flag)
    K2x = h * K2x
    K2y = h * K2y

    K3x, K3y = dx(x0 + K2x/2, y0 + K2y/2, flag)
    K3x = h * K3x
    K3y = h * K3y

    K4x, K4y = dx(x0 + K3x, y0 + K3y, flag)
    K4x = h * K4x
    K4y = h * K4y

    Kx = (1 / 6) * (K1x + 2 * K2x + 2 * K3x + K4x)
    x0 = x0 + Kx

    Ky = (1 / 6) * (K1y + 2 * K2y + 2 * K3y + K4y)
    y0 = y0 + Ky

    return x0, y0

def rk4_loop(x0, y0, h, n):
    '''
        Calculate the approximate value of CO2-concentration of the greenhouse air
        and CO2-concentration of the top greenhouse air at t + h, t + 2h, ..., t + n*h

        :param x0: CO2-concentration of the greenhouse air at time start t
        :param y0: CO2-concentration of the top greenhouse air at time start t
        :param h: step approximation
        :param n: the number of times you want to get data (means you can still get data after
                  n*h (minutes))
    '''

    print("------RUNGE-KUTTA 4 RESULT------")

    file_meteo = open("meteo.csv")
    file_meteo.readline()
    file_ghc = open("Greenhouse_climate.csv")
    file_ghc.readline()
    file_vip = open("vip.csv")
    file_vip.readline()

    CO2_air = x0
    CO2_top = y0
    step_next = 1
    step = 0.1
    n *= 300
    i = 0
    t0 = 0
    getData = False
    CO2_air_real_data = 427

    while i < n:
        if i // 300 >= step_next:
            step_next += 1

            line_meteo = file_meteo.readline()
            line_ghc = file_ghc.readline()
            line_vip = file_vip.readline()
            elements_vip = line_vip.split(',')

            while "NaN" in line_meteo or "NaN" in line_ghc or elements_vip[1] == "NaN":
                line_meteo = file_meteo.readline()
                line_ghc = file_ghc.readline()
                line_vip = file_vip.readline()
                elements_vip = line_vip.split(',')

            elements_vip = line_vip.split(',')
            elements_meteo = line_meteo.split(',')
            elements_ghc = line_ghc.split(',')

            init.v_wind = float(elements_meteo[10])
            init.T_out = float(elements_meteo[8])
            init.T_air = float(elements_ghc[9])
            init.U_thscr = float(elements_ghc[3]) / 100
            init.U_roof = (float(elements_ghc[10]) + float(elements_ghc[11])) / 200
            init.T_can = init.T_air + 15
            init.T_top = init.T_air + 1

            if elements_vip[0] == "NaN":
                init.U_extCO2 = np.random.uniform(0, 0.3)
            elif init.convert_mgm3_to_ppm(CO2_air) < float(elements_vip[0]):
                init.U_extCO2 = 1
            else:
                init.U_extCO2 = 0

            CO2_air_real_data = float(elements_ghc[2])

            getData = True

        if getData:
            randomU()
            getData = False

        CO2_air_step, CO2_top_step = rk4(dx, CO2_air, CO2_top, step, getData)
        CO2_air = CO2_air_step
        CO2_top = CO2_top_step

        if t0 % (h * 60) == 0:

            print("CO2_real: ", CO2_air_real_data)
            list_real_CO2_air_data.append(CO2_air_real_data)

            err = abs(CO2_air_real_data - CO2_air)
            list_err_rk4.append(init.convert_mgm3_to_ppm(err))

            list_CO2_air_rk4.append(init.convert_mgm3_to_ppm(CO2_air))
            list_CO2_top_rk4.append(init.convert_mgm3_to_ppm(CO2_top))

            print("CO2_air at t +", t0 // 60, "=", init.convert_mgm3_to_ppm(CO2_air))
            print("CO2_air at t +", t0 // 60, "=", init.convert_mgm3_to_ppm(CO2_top))
            print("------------------")

        i += 1
        t0 += 1

    file_meteo.close()
    file_ghc.close()
    file_vip.close()


CO2_air_0 = init.convert_ppm_to_mgm3(427)
CO2_top_0 = init.convert_ppm_to_mgm3(427)

# CO2_air_diff, CO2_top_diff = dx(CO2_air_0, CO2_top_0)
# print("CO2_air_diff =", init.convert_mgm3_to_ppm(CO2_air_diff))
# print("CO2_top_diff =", init.convert_mgm3_to_ppm(CO2_top_diff))

step = 5
n = 1000

# euler_loop(CO2_air_0, CO2_top_0, step, n)
rk4_loop(CO2_air_0, CO2_top_0, step, n)

err_mean = 0
for i in range(len(list_err_rk4)):
    err_mean += list_err_rk4[i]

err_mean = err_mean / n
print("Err Mean= ", err_mean)

x_real = [i for i in range(n)]

plt.figure()
plt.title("Runge-Kutta Result")

plt.subplot(1, 2, 1)
plt.plot(x_real, list_real_CO2_air_data, label='CO2_real_data')
plt.plot(x_real, list_CO2_air_rk4, label='CO2_air')
plt.plot(x_real, list_CO2_top_rk4, label='CO2_top')
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(x_real, list_err_rk4, label='Error')
plt.legend()

# plt.title("Euler Result")
# plt.subplot(1, 2, 1)
# plt.plot(x_real, list_real_CO2_air_data, label='CO2_real_data')
# plt.plot(x_real, list_CO2_air_euler, label='CO2_air')
# plt.plot(x_real, list_CO2_top_euler, label='CO2_top')
# plt.legend()
#
# plt.subplot(1, 2, 2)
# plt.plot(x_real, list_err_euler, label='Error')
# plt.legend()

plt.xlabel('Time')
plt.ylabel('CO2-concentration (ppm)')

plt.show()
