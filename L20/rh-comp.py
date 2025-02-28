from lagrange_gen import Experiment
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import os

import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredText

'''
 calc_a function and associated constants adapted from  https://github.com/icepack/icepack/blob/master/src/icepack/models/viscosity.py
'''

MAX_TEMP = 0
MIN_TEMP = -30

YEAR = 365.25 * 24 * 60 * 60 #seconds in a year

CRITICAL_TEMP  = 263.15
IDEAL_GAS = 8.3144621e-3

A0_cold = 3.985e-13 * YEAR * 1.0e18  # 1 / (MPa^3 yr)
A0_warm = 1.916e3 * YEAR * 1.0e18
Q_cold = 60  # kJ / mol
Q_warm = 139

strain_over_time = 10e-5  # per year

Eij_factor = 1

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


def glen_law(temp, t, n = 3):
    '''
    Temp in celsius
    '''
    e = calc_a(c_to_K(temp)) * (t ** n)
    return e



# https://github.com/nicholasmr/specfab/blob/main/src/deformationmodes.f90

'''
diag([b**((1+r)/2), b**((1-r)/2), b**(-1)])
'''

def ss_strain(T,angle):
    return T **( (np.deg2rad(angle) + 1) / 2)



def plot_experiment_enhancement(ex, T):
    '''
    ex is a list of Experiment objects which have already been generated as an .nc file
    T is the time the strain was over

    '''



    dpi, scale = 200, 3.3
    fig = plt.figure(figsize=(4*scale,3*scale))
    gs = gridspec.GridSpec(2,2, height_ratios=[1,1], width_ratios=[1,1])
    #gs.update(left=-0.03, right=1-0.06/3, top=0.97, bottom=0.20, wspace=0.015*18, hspace=0.35)


    ax_Elin    = plt.subplot(gs[0, 0])
    ax_Enlin   = plt.subplot(gs[0, 1])
    ax_Elin_E    = plt.subplot(gs[1, 0])
    ax_Enlin_E   = plt.subplot(gs[1, 1])
  
    lw0, lw1, lw2 = 2.5,2.25,2.0


    colorbar_made = False
    current_slopes = {
        ax_Elin: [],
        ax_Enlin: [],
        ax_Elin_E: [],
        ax_Enlin_E: []
    }
    ax_Elin_data =  None
    ax_Enlin_data = None
    ax_Elin_E_data   =  None
    ax_Enlin_E_data = None

    for x in ex:
        
        ax_Elin_data = plot_enhancements(ax_Enlin, 'nonlinear_enhancement', x, T, current_slopes, use_Eij = False)
        ax_Enlin_data = plot_enhancements(ax_Elin,  'linear_enhancement', x, T, current_slopes, use_Eij = False)
        ax_Elin_E_data = plot_enhancements(ax_Enlin_E, 'nonlinear_enhancement', x, T, current_slopes)
        ax_Enlin_E_data = plot_enhancements(ax_Elin_E,  'linear_enhancement', x, T, current_slopes)



    ax_Enlin.set_title(r'Nonlinear Sachs')
    add_slope(ax_Enlin, current_slopes)
    fig.colorbar(ax_Elin_data, label = "Temperature (°C)")

    ax_Elin.set_title(r'Linear mixed Taylor--Sachs')
    add_slope(ax_Elin, current_slopes)
    fig.colorbar(ax_Enlin_data, label = "Temperature (°C)")

    ax_Enlin_E.set_title(r'Nonlinear Sachs with E')
    add_slope(ax_Enlin_E, current_slopes)
    fig.colorbar(ax_Elin_E_data, label = "Temperature (°C)")

    ax_Elin_E.set_title(r'Linear mixed Taylor--Sachs with E')
    add_slope(ax_Elin_E, current_slopes)
    fig.colorbar(ax_Enlin_E_data, label = "Temperature (°C)")



    os.makedirs("output", exist_ok=True)
    fout = f'output/cache.png'
    print('Saving %s'%(fout))
    plt.savefig(fout, dpi=dpi)
    plt.close('all')




def add_slope(ax, slope_dict):
    ax.text(.05,.95, "Average n = " + str(round(np.mean(slope_dict[ax]), 2)) + 
            "\n Min n = " + str(round(np.min(slope_dict[ax]), 2)) +
            "\n Max n = " + str(round(np.max(slope_dict[ax]), 2)) , 
            
            transform=ax.transAxes, horizontalalignment='left', verticalalignment='top',
            bbox=dict(facecolor='white', alpha=0.5, edgecolor='black', boxstyle='round,pad=0.25'))




def plot_enhancements(ax, Eij, ex, T, slopes, use_Eij = True):
    df = ex.get_dataframe()
    df.strain = np.abs(df.strain)
    print(df[Eij])
    if ex.exptype != "ss":
        df['tau'] = df.strain * T
    else:
        df['tau'] = ss_strain(T, df.strain)
        
    #df['nondistinct_tau'] = np.array([df.tau[len(df.tau) - 1]] * len(df.strain))


    x_ = df['tau']
    if use_Eij == "exponent":
        ax.set_ylabel('log(Glens where n=3+E)')
        df['glens'] = glen_law(ex.temp, df['tau'], n=4 * df[Eij])
        y_ =  df['glens']
    elif use_Eij:
        ax.set_ylabel('log(Glens * E)')
        df['glens'] = glen_law(ex.temp, df['tau'])
        y_ =  df['glens'] * (df[Eij] * Eij_factor)
    else:
        ax.set_ylabel('log(Glens)')
        df['glens'] = glen_law(ex.temp, df['tau'])
        y_ =  df['glens'] #* Eij * Eij_factor

    
    data = ax.scatter(x_, y_, label=str(ex.temp) + "°C", c=[ex.temp]* len(df.strain), s = 2, norm=colors.Normalize(MIN_TEMP, MAX_TEMP))

    if ex.exptype == "ss":
        ax.set_xlabel('Target Angle')
    else:
        ax.set_xlabel('Target Change')


    ax.set_xlabel("log(Tau)")

    ax.grid()
    #ax.legend(fontsize=7)
    ax.set_xscale('log')
    ax.set_yscale('log')

    m, b = np.polyfit(np.log(x_[1:]), np.log(y_[1:]), 1)
    slopes[ax].append(m)

    return data



scope = []

EXP = ["uc", "ss", "cc", "ue"]
CC = ["cc"]
SS = ["ss"]
UE = ["ue"]
UC = ["uc"]

TEMPS = []
for x in range(MIN_TEMP, MAX_TEMP, 2):
    TEMPS.append(x)


for x in SS:
    for tem in TEMPS:
        print(x, tem)
        
        if x != "ss" :
            e1 = Experiment(x, "zz", temp = tem) 
        else:
            e1 = Experiment(x, "xz", temp = tem) 

        print(f"GAMMA: {e1.Gamma}, LAMD: 0.15, TEMP: {e1.temp}, EXP: {e1.exptype}")
        scope.append(e1)


#e1 = Experiment("ss", "xz", temp = -30) 

plot_experiment_enhancement(scope, strain_over_time)