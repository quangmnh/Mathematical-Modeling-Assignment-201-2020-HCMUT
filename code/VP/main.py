import init
import random
import os
import sys
import math
import random
import numpy as np
from matplotlib import pyplot as plt

sys.stdout = open('logger.txt','w')
''' ------------------------ '''
''' System params            '''
''' ------------------------ '''
step = 0.1      #Bước xấp xỉ
t=500           #Thời gian chạy model
VP_air_0 = init.VP_sat(19.8999999966472)*0.804
VP_top_0 = init.VP_sat(18.8999999966472)*0.804


def f_func():
    return (init.MV_can_air() + \
           init.MV_pad_air() + \
           init.MV_fog_air() + \
           init.MV_blow_air()- \
           init.MV_air_thscr() - \
           init.MV_air_top() - \
           init.MV_air_out() - \
           init.MV_air_out_pad() - \
           init.MV_air_mech())/init.cap_VP_air()


def g_func():
    return (init.MV_air_top() - \
           init.MV_top_cov_in() - \
           init.MV_top_out())/init.cap_VP_top()

air_top1 = []
top_cov = []
top_out =[]
can_air = []
air_top = []
air_out = []
air_thscr = []

def control_val():
    valve = 0.5
    init.U_fog = 0.5
    init.U_mech_cool = 0.5
    init.U_pad = 0.5
    init.U_blow = 0.2
    init.U_side = 0.5
    init.U_vent_forced = 0.1
    init.U_ThScr = 0.5
    init.U_shscr = 0.5
    init.U_thscr = init.U_ThScr


def param_update(elements,elements2):
    init.v_wind = float(elements[10])
    init.T_out = float(elements[8])
    init.R_can = 0.65 *float(elements[2])
    init.T_air = float(elements2[9])
    init.U_ThScr = float(elements2[3])/100
    init.U_thscr = float(elements2[3])/100
    init.U_roof = abs((float(elements2[10])+float(elements2[11])))/200
    init.CO2_air = float(elements2[2])
    init.T_top = init.T_air - 1
    init.T_thscr = init.T_air + 1
    init.T_mech_cool = init.T_air - 1
    init.T_top_cov_in = init.T_top - 1
    init.T_can = init.T_air+15
    init.T_cov_in = init.T_top-1
    init.T_mean_air =0.5*(init.T_air+init.T_top)
    init.VP_out = float(elements[7])*init.VP_Out()/100
def dx(x0, y0):
    '''
    Calculate the rate of change in CO2 concentration in greenhouse and top of greenhouse

    :param x0: Vapour Pressure of the greenhouse air at time start t
    :param y0: Vapour Pressure of the top greenhouse air at time start t
    '''
    init.VP_air = x0
    init.VP_top = y0
    VP_air = f_func()
    VP_top = g_func()

    return VP_air, VP_top
def euler(x0, y0, h):
    '''
    Calculate the approximate value of Vapour Pressure of the greenhouse air
    and Vapour Pressure of the top greenhouse air

    :param dx: the callable function of the ODE system
    :param x0: Vapour Pressure of the greenhouse air at time start t
    :param y0: Vapour Pressure of the top greenhouse air at time start t
    :param h: step approximation
    '''
    Kx, Ky = dx(x0, y0)
    Kx = h * Kx
    Ky = h * Ky

    x0 = x0 + Kx
    y0 = y0 + Ky

    return x0, y0
def euler_solve(x0, y0, h):
    file = open("meteo.csv")                
    file.readline()
    n=h//5  #Tính số lần xấp xỉ để đủ 5 phút
    file2 = open("Greenhouse_climate.csv")
    file2.readline()
    for i in range(n):
        line = file.readline()
        while "NaN" in line:
            line = file.readline()
        elements = line.split(',')


        line2 = file2.readline()
        while "NaN" in line2:
            line2 = file2.readline()
        elements2 = line2.split(',')

        #Update các giá trị mỗi 5 phút
        param_update(elements,elements2)
        
        #Xấp xỉ cho 1 khoảng thời gian 5 phút
        for j in range(int(300//step)):
            x0,y0=euler(x0,y0,step)

        #Giới hạn giá trị RH lại nếu lớn hơn 100%
        if x0>init.VP_Air():
            x0 = init.VP_Air()
        if y0>init.VP_Top():
            y0 = init.VP_Top()

        print("--------------------------------")
        print("VP Air at t +", i * 5, "= ", x0)
        print("VP Top at t +", i * 5, "= ", y0)
        print("Relative humidity Air at t +", i * 5, "= ", x0/init.VP_Air()*100)
        print("relative humidity Top at t +", i * 5, "= ", y0/init.VP_Top()*100)
        
        eu_air.append(x0/init.VP_Air()*100)
        eu_top.append(y0/init.VP_Top()*100)

    file.close()
    file2.close()
    return x0, y0
def rk4(x0,y0,h):
    '''
    Calculate the approximate value of CO2-concentration of the greenhouse air
    and CO2-concentration of the top greenhouse air

    :param dx: the callable function of the ODE system
    :param x0: CO2-concentration of the greenhouse air at time start t
    :param y0: CO2-concentration of the top greenhouse air at time start t
    :param h: step approximation
    :param flag: To be used for printing result to log file and saving to data for plotting
    '''

    K1x, K1y = dx(x0, y0)
    K1x = h * K1x
    K1y = h * K1y

    K2x, K2y = dx(x0 + K1x/2, y0 + K1y/2)
    K2x = h * K2x
    K2y = h * K2y

    K3x, K3y = dx(x0 + K2x/2, y0 + K2y/2)
    K3x = h * K3x
    K3y = h * K3y

    K4x, K4y = dx(x0 + K3x, y0 + K3y)
    K4x = h * K4x
    K4y = h * K4y

    Kx = (1 / 6) * (K1x + 2 * K2x + 2 * K3x + K4x)
    x0 = x0 + Kx

    Ky = (1 / 6) * (K1y + 2 * K2y + 2 * K3y + K4y)
    y0 = y0 + Ky

    return x0, y0
def rk4_solve(x0, y0, h):
    file = open("meteo.csv")
    file.readline()
    n=h//5 #Tính số lần xấp xỉ để đủ 5 phút
    file2 = open("Greenhouse_climate.csv")
    file2.readline()
    for i in range(n):

        line = file.readline()
        while "NaN" in line:
            line = file.readline()
        elements = line.split(',')

        line2 = file2.readline()
        while "NaN" in line2:
            line2 = file2.readline()
        elements2 = line2.split(',')
        
        #Update các giá trị mỗi 5 phút
        param_update(elements,elements2)
        
        #Xấp xỉ cho mỗi khoảng thời gian 5 phút
        for j in range(int(300//step)):
            x0,y0=rk4(x0,y0,step)
        
        #Giới hạn giá trị RH lại nếu lớn hơn 100%
        if x0>init.VP_Air():
            x0 = init.VP_Air()
        if y0>init.VP_Top():
            y0 = init.VP_Top()

        print("--------------------------------")
        print("VP Air at t +", i * 5, "= ", x0)
        print("VP Top at t +", i * 5, "= ", y0)
        print("Relative Humidity Air at t +", i * 5, "= ", x0/init.VP_Air()*100)
        print("Relative Humidity Top at t +", i * 5, "= ", y0/init.VP_Top()*100)
        rk_air.append(x0/init.VP_Air()*100)
        rk_top.append(y0/init.VP_Top()*100)
    file.close()
    return x0, y0

def cons(x):
    if x>1:
        return 1
    if x<0:
        return 0
    return x

def random_val():
    init.U_fog = cons(random.gauss(0.6,0.2))
    init.U_mech_cool = cons(random.gauss(0.6,0.2))
    init.U_pad = cons(random.gauss(0.6,0.2))
    init.U_blow = cons(random.gauss(0.6,0.2))
    init.U_side = cons(random.gauss(0.6,0.2))
    init.U_vent_forced = cons(random.gauss(0.6,0.2))
    init.U_ThScr = cons(random.gauss(0.6,0.2))
    init.U_shscr = cons(random.gauss(0.6,0.2))
    init.U_thscr = init.U_ThScr

control_val()

eu_air =[]
rk_air = []

eu_top = []
rk_top = []

eu_air_err = []
rk_air_err = []

print(" EULER RESULT: ")
euler_solve(VP_air_0, VP_top_0, t)

print(" RUNGE-KUTTA 4 RESULT: ")
rk4_solve(VP_air_0, VP_top_0, t)
time = []
for i in range(t//5):
    time.append(i*5)
real =[]
file = open("Greenhouse_climate.csv")
file.readline()

for i in range(t//5):
    line = file.readline()

    while "NaN" in line:
        line = file.readline()

    elements = line.split(',')
    real.append(float(elements[8]))
print(len(real))
print(len(eu_air))
for i in range(t//5):
    eu_air_err.append(abs(real[i]-eu_air[i]))
    rk_air_err.append(abs(real[i]-eu_top[i]))
print(sum(eu_air_err)/len(eu_air_err))
print(sum(rk_air_err)/len(rk_air_err))

plt.figure()
plt.subplot(1,2,1)
plt.title("Euler results")
plt.plot(time,real,label ='Observed data')
plt.plot(time,eu_air,label ='Vapour pressure Air')
plt.plot(time,eu_top,label = 'Vapour pressure Top')
plt.legend()
plt.subplot(1,2,2)
plt.title("Error of Euler method")
plt.plot(time,eu_air_err,label ='Error')
plt.legend()

plt.figure()
plt.subplot(1,2,1)
plt.title("Runge-kutta results")
plt.plot(time,real,label ='Observed data')
plt.plot(time,rk_air,label = 'Vapour pressure Air')
plt.plot(time,rk_top,label = 'Vapour pressure Top')
plt.legend()
plt.subplot(1,2,2)
plt.title("Error of Runge-kutta method")
plt.plot(time,rk_air_err,label ='Error')
plt.legend()


plt.legend()
plt.show()

