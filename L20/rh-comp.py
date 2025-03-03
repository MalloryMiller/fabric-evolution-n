
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import os
from utils import *

import matplotlib.gridspec as gridspec

'''
 calc_a function and associated constants adapted from  https://github.com/icepack/icepack/blob/master/src/icepack/models/viscosity.py
'''



def plot_all_enhancement(ex, t, cutoff = True):
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
    
    current_slopes = {}

    cutoff_title = ""
    if cutoff:
        cutoff_title = ""
  
    plot_enhancement(ax_Elin, fig, "nonlinear_enhancement", ex, t, current_slopes, "Nonlinear Sachs", False, cutoff = cutoff)
    plot_enhancement(ax_Enlin, fig, "linear_enhancement", ex, t, current_slopes, "Linear mixed Taylor--Sachs", False, cutoff = cutoff)
    plot_enhancement(ax_Elin_E, fig, "nonlinear_enhancement", ex, t, current_slopes, "Nonlinear Sachs", cutoff = cutoff)
    plot_enhancement(ax_Enlin_E, fig, "linear_enhancement", ex, t, current_slopes, "Linear mixed Taylor--Sachs", cutoff = cutoff)


    os.makedirs("output", exist_ok=True)
    fout = f'output/cache.png'
    print('Saving %s'%(fout))
    plt.savefig(fout, dpi=dpi)
    plt.close('all')



def plot_enhancement(ax, fig, Ei, ex, t, slopes, title, use_Eij = True, cutoff = True):
    slopes[ax] = []

    for x in ex:
        ax_data = generate_plot(ax, Ei, x, t, slopes, use_Eij = use_Eij, cutoff=cutoff)

    ax.set_title(title)
    add_slope(ax, slopes)
    fig.colorbar(ax_data, label = "Temperature (°C)")





def add_slope(ax, slope_dict):
    ax.text(.05,.95, "Average n = " + str(round(np.mean(slope_dict[ax]), 2)) + 
            "\n Min n = " + str(round(np.min(slope_dict[ax]), 2)) +
            "\n Max n = " + str(round(np.max(slope_dict[ax]), 2)), 
            
            transform=ax.transAxes, horizontalalignment='left', verticalalignment='top',
            bbox=dict(facecolor='white', alpha=0.5, edgecolor='black', boxstyle='round,pad=0.25'))




def generate_plot(ax, Eij, ex, t, slopes, use_Eij = True, cutoff = True):
    df = ex.get_dataframe()
    df.strain = np.abs(df.strain)

    y_ = df.strain / t # strain rate 
    if use_Eij:
        ax.set_ylabel('log(Strain Rate * E)')
        x_ = reverse_glens(df.strain, ex.temp, t, E=df[Eij])
    else:
        ax.set_ylabel('log(Strain Rate)')
        x_ = reverse_glens(df.strain, ex.temp, t)

    temps = [ex.temp] * len(df.strain)


    if cutoff:
        x_ = cut_to_range(x_, df[Eij])
        y_ = cut_to_range(y_, df[Eij])
        temps = cut_to_range(temps, df[Eij])

        
    
    data = ax.scatter(x_, y_, label=str(ex.temp) + "°C", c=temps, s = 2, norm=colors.Normalize(MIN_TEMP, MAX_TEMP))

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



#e1 = Experiment("ss", "xz", temp = -30) 

plot_all_enhancement(scope, strain_over_time)

