import numpy as np
import pandas as pd

from lagrange_gen import Experiment

MAX_TEMP = 0
MIN_TEMP = -30

YEAR = 365.25 * 24 * 60 * 60 #seconds in a year

CRITICAL_TEMP  = 263.15
IDEAL_GAS = 8.3144621e-3

A0_cold = 3.985e-13 * YEAR * 1.0e18  # 1 / (MPa^3 yr)
A0_warm = 1.916e3 * YEAR * 1.0e18
Q_cold = 60  # kJ / mol
Q_warm = 139

strain_over_time = 2e-1  # per year

Eij_factor = 1

Eij_cutoff = 0.00001


scope = []

EXP = ["uc", "ss", "cc", "ue"]
CC = ["cc"]
SS = ["ss"]
UE = ["ue"]
UC = ["uc"]

TEMPS = []
for x in range(MIN_TEMP, MAX_TEMP+1, 2):
    TEMPS.append(x)


for x in CC:
    for tem in TEMPS:
        
        if x != "ss" :
            e1 = Experiment(x, "zz", temp = tem) 
        else:
            e1 = Experiment(x, "xz", temp = tem) 

        print(f"GAMMA: {e1.Gamma}, LAMD: 0.15, TEMP: {e1.temp}, EXP: {e1.exptype}")
        scope.append(e1)





def get_range(Es):
    i = 1
    to_test = []
    while i < len(Es):
        to_test.append(Es[i-1] - Es[i])
        if np.abs(Es[i-1] - Es[i]) < Eij_cutoff:
            return i
        i += 1

    print("NEVER STOPPED: ", np.min(to_test))
    return len(Es)-1


def cut_to_range(Cut, Es):
    '''
    Cut : the list being cut down based on Es
    Es  : where the range being cut of is derived
    '''
    i = get_range(Es)
    return Cut[:i]



def get_balanced_average(x, y, temp = 0):
    x = np.sort(x)
    y = np.sort(y)
    x_range = x[-1] - x[0]
    average = 0
    
    i = 1
    while i < len(x):
        rise = y[i] - y[i-1]
        run  = x[i] - x[i-1]

        average += (rise / run) * (run / x_range)

        i += 1



    results = {
        "ns": [average], 
        "temp": [temp], 
        "weight": [1], 
        }

    return pd.DataFrame(

        data = results
        
    )


def get_individual_slopes(x, y, temp):
    x = np.copy(x)
    y = np.copy(y)
    temp = temp[:len(x) - 1]
    slopes = []
    x_range = x[-1] - x[0]
    weights = []

    
    i = 1
    while i < len(x):
        rise = y[i] - y[i-1]
        run  = x[i] - x[i-1]

        slopes.append(rise / run)
        weights.append(run/x_range)


        i += 1


    results = pd.DataFrame(

        data = {

            "ns": slopes, 
            "temp": temp, 
            "weight": weights, 

            }
        
    )

    return results





def calc_a(temp):
    '''
    Temp in Kelvin
    '''
    A0 = A0_warm
    Q = Q_warm
    if temp < CRITICAL_TEMP:
        A0 = A0_cold
        Q = Q_cold

    return A0 * np.exp(-Q / (IDEAL_GAS * temp))





def K_to_c(K):
    '''
    Kelvin to celsius
    '''
    return K - 273.15

def c_to_K(c):
    '''
    Celsius to Kelvin
    '''
    return c + 273.15


def glen_law(temp, tau, n = 3):
    '''
    Temp in celsius
    '''
    e = calc_a(c_to_K(temp)) * (tau ** n)
    return e



def reverse_glens(strain, temp, time, E=1, n=3):
    '''
    Temp in celsius
    '''
    return (strain / (calc_a(c_to_K(temp)) * time * E)) ** (1/n)




# https://github.com/nicholasmr/specfab/blob/main/src/deformationmodes.f90

'''
diag([b**((1+r)/2), b**((1-r)/2), b**(-1)])
'''

def ss_strain(T,angle):
    return (T **( (np.deg2rad(angle) + 1) / 2))
