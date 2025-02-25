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

strain_over_time = 10e-5

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
    return K - 273.15

def c_to_K(c):
    return c + 273.15


def glen_law(temp, t, n = 3):
    '''
    Temp in celsius
    '''
    e = calc_a(c_to_K(temp)) * (t ** n)
    return e



def plot_experiment_enhancement(ex, T):
    '''
    ex is an Experiment object which already has generated an .nc file

    '''



    dpi, scale = 200, 3.3
    fig = plt.figure(figsize=(4/2*2.1*scale,2.3*scale/1.7))
    gs = gridspec.GridSpec(1,2, height_ratios=[.5], width_ratios=[1,1])
    #gs.update(left=-0.03, right=1-0.06/3, top=0.97, bottom=0.20, wspace=0.015*18, hspace=0.35)

    ax_Elin    = plt.subplot(gs[0, 0])
    ax_Enlin   = plt.subplot(gs[0, 1])
  
    lw0, lw1, lw2 = 2.5,2.25,2.0


    colorbar_made = False

    for x in ex:

        df = x.get_dataframe()

        def plot_enhancements(ax, Eij, df, colorbar_made = colorbar_made):
            if x.exptype == "ss":
                df.strain = df.strain / 90
            else:
                df.strain = np.abs(df.strain)
            print(df.strain)

            df['tau'] = df.strain / T
            df['nondistinct_tau'] = np.array([df.tau[len(df.tau) - 1]] * len(df.strain))

            df['glens'] = glen_law(x.temp, df['tau'])

            x_ = df['glens'] *Eij
            y_ = df['tau']

            
            #ax.semilogy(steps, Eij[:,0], '-', c=sfplt.c_red,   label=lblm(0,0), lw=lw0)
            data = ax.scatter(x_, y_, label=str(x.temp) + "°C", c=[x.temp]* len(df.strain), s = 2, norm=colors.Normalize(MIN_TEMP, MAX_TEMP))

            if x.exptype == "ss":
                ax.set_xlabel('Target Angle')
            else:
                ax.set_xlabel('Target Change')

            ax.set_ylabel("Tau")

            ax.set_xlabel('GlensxE')
            ax.grid()
            #ax.legend(fontsize=7)
            ax.set_xscale('log')
            ax.set_yscale('log')
            if not colorbar_made:
                fig.colorbar(data, label = "Temperature (C)")
                colorbar_made = True


            return colorbar_made
        

        colorbar_made = plot_enhancements(ax_Enlin, df.nonlinear_enhancement, df)
        ax_Enlin.set_title(r'Nonlinear Sachs')

        plot_enhancements(ax_Elin,  df.linear_enhancement, df)
        ax_Elin.set_title(r'Linear mixed Taylor--Sachs')



    os.makedirs("output", exist_ok=True)
    fout = 'output/cache.png'
    print('Saving %s'%(fout))
    plt.savefig(fout, dpi=dpi)
    plt.close('all')


scope = []

EXP = "cc"
TEMPS = []
for x in range(MIN_TEMP, MAX_TEMP, 2):
    TEMPS.append(x)



for tem in TEMPS:
    print(EXP, tem)
    
    if EXP != "ss" :
        e1 = Experiment(EXP, "zz", temp = tem) 
    else:
        e1 = Experiment(EXP, "xz", temp = tem) 

    print(f"GAMMA: {e1.Gamma}, LAMD: 0.15, TEMP: {e1.temp}, EXP: {e1.exptype}")
    scope.append(e1)


#e1 = Experiment("ss", "xz", temp = -30) 

plot_experiment_enhancement(scope, strain_over_time)