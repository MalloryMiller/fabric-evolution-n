
from matplotlib.cm import ScalarMappable
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mp
import numpy as np
import os
from utils import *

import matplotlib.gridspec as gridspec


CMAP_TEMP = mp.cm.plasma
CMAP_STRAIN_RATE = mp.cm.copper

'''
 calc_a function and associated constants adapted from  https://github.com/icepack/icepack/blob/master/src/icepack/models/viscosity.py
'''



def plot_all_enhancement(ex, t, cutoff = True, balanced_average = True, fname = "cache"):
    '''
    ex is a list of Experiment objects which have already been generated as an .nc file
    T is the time the strain was over

    '''

    dpi, scale = 200, 3.3
    fig = plt.figure(figsize=(4*scale,8.5*scale))
    gs = gridspec.GridSpec(5,2, height_ratios=[1,1,.75,.75,.75], width_ratios=[1,1])
    #gs.update(left=-0.03, right=1-0.06/3, top=0.97, bottom=0.20, wspace=0.015*18, hspace=0.35)


    ax_Elin               = plt.subplot(gs[0, 0])
    ax_Enlin              = plt.subplot(gs[0, 1])
    ax_Elin_E             = plt.subplot(gs[1, 0])
    ax_Enlin_E            = plt.subplot(gs[1, 1])
    ax_Elin_n_temp        = plt.subplot(gs[2, 0])
    ax_Enlin_n_temp       = plt.subplot(gs[2, 1])
    ax_Elin_n_rate        = plt.subplot(gs[3, 0])
    ax_Enlin_n_rate        = plt.subplot(gs[3, 1])
    ax_Elin_n_temp_extr    = plt.subplot(gs[4, 0])
    ax_Enlin_n_temp_extr   = plt.subplot(gs[4, 1])
    
    current_slopes = {}

    cutoff_title = ""
    if cutoff:
        cutoff_title = ""
  
    plot_enhancement(ax_Elin, fig, "nonlinear_enhancement", ex, t, current_slopes, "Nonlinear Sachs", 
                     False, cutoff = cutoff, balanced_average = balanced_average)
    plot_enhancement(ax_Enlin, fig, "linear_enhancement", ex, t, current_slopes, "Linear mixed Taylor--Sachs", 
                     False, cutoff = cutoff, balanced_average = balanced_average)
    plot_enhancement(ax_Elin_E, fig, "nonlinear_enhancement", ex, t, current_slopes, "Nonlinear Sachs", 
                     cutoff = cutoff, balanced_average = balanced_average)
    plot_enhancement(ax_Enlin_E, fig, "linear_enhancement", ex, t, current_slopes, "Linear mixed Taylor--Sachs", 
                     cutoff = cutoff, balanced_average = balanced_average)


    n_hist(fig, ax_Elin_n_temp, current_slopes[ax_Elin_E],"Nonlinear Sachs")
    n_hist(fig, ax_Enlin_n_temp, current_slopes[ax_Enlin_E], "Linear mixed Taylor--Sachs")

    n_hist(fig, ax_Elin_n_rate, current_slopes[ax_Elin_E],"Nonlinear Sachs", which_ns = 'avg_st_rate')
    n_hist(fig, ax_Enlin_n_rate, current_slopes[ax_Enlin_E], "Linear mixed Taylor--Sachs", which_ns = 'avg_st_rate')

    extremes_hist(fig, ax_Elin_n_temp_extr, current_slopes[ax_Elin_E],"Nonlinear Sachs")
    extremes_hist(fig, ax_Enlin_n_temp_extr, current_slopes[ax_Enlin_E], "Linear mixed Taylor--Sachs")

    os.makedirs("output", exist_ok=True)
    fout = f'output/{fname}.png'
    print('Saving %s'%(fout))
    plt.savefig(fout, dpi=dpi)
    plt.close('all')



def plot_enhancement(ax, fig, Ei, ex, t, slopes, title, use_Eij = True, cutoff = True, balanced_average = True):
    results = {
        "ns": [], 
        "temp": [], 
        "weight": []
        }

    slopes[ax] = pd.DataFrame(

        data = results
        
    )

    for x in ex:
        ax_data = generate_plot(ax, Ei, x, t, slopes, use_Eij = use_Eij, cutoff=cutoff, balanced_average = balanced_average)

    ax.set_title(title)
    add_slope(ax, slopes[ax])
    fig.colorbar(ax_data, label = "Temperature (°C)")





def add_slope(ax, slope_dict):
    ax.text(.05,.95, "Average n = " + str(round(np.mean(slope_dict.ns), 2)) + 
            "\n Min n = " + str(round(np.min(slope_dict.ns), 2)) +
            "\n Max n = " + str(round(np.max(slope_dict.ns), 2)), 
            
            transform=ax.transAxes, horizontalalignment='left', verticalalignment='top',
            bbox=dict(facecolor='white', alpha=0.5, edgecolor='black', boxstyle='round,pad=0.25'))
    

def extremes_hist(fig, ax, ns, Etype, which_ns = 'temp'):
    datas = []

    min_ns = np.min(ns[which_ns].unique())
    max_ns = np.max(ns[which_ns].unique())
    temps = [min_ns, max_ns]
    min_datas = np.min([np.min(ns.ns[ns[which_ns] == min_ns]), np.min(ns.ns[ns[which_ns] == max_ns])])
    max_datas = np.max([np.max(ns.ns[ns[which_ns] == min_ns]), np.max(ns.ns[ns[which_ns] == max_ns])])
    datas.append(ns.ns[ns[which_ns] == min_ns])
    datas.append(ns.ns[ns[which_ns] == max_ns])
    bin_count = 20

    bins = []

    for x in range(bin_count + 1):
        bins.append(min_datas + (x * ((max_datas - min_datas) / bin_count)))

    print("TEMPS", temps)


    ax.hist(datas[0], density=True, rwidth=0.95, stacked = False, alpha=.75, color="darkslateblue", label=str(temps[0]) + " °C", bins = bins)
    ax.hist(datas[1], density=True, rwidth=0.95, stacked = False, alpha=.75, color="gold", label=str(temps[1]) +  " °C", bins = bins)

    #ax.bar(n, bins, histtype='bar', rwidth=0.95, stacked = True)
    title = str(Etype) + " n Values"
    ax.set_title(title)
    ax.set_xlabel("n Value")
    ax.set_ylabel("Count")
    ax.legend()






def n_hist(fig, ax, ns, Etype, which_ns = 'temp'):
    datas = []
    temps = []

    for t in ns[which_ns].unique():
        temps.append(t)
        datas.append(ns.ns[ns[which_ns] == t])
        
    min_ns = np.min(temps)
    max_ns = np.max(temps)

    CMAP = CMAP_TEMP
    norm = colors.Normalize(min_ns, max_ns)
    c_label = "Temperature (°C)"
    if which_ns == 'avg_st_rate':
        CMAP = CMAP_STRAIN_RATE
        norm = colors.LogNorm(min_ns, max_ns)
        c_label = "log(Strain Rate x E)"


    h, bins, patches = ax.hist(datas, histtype='bar', rwidth=0.95, stacked = True)
    for i, patch in enumerate(patches):
        for bar in patch:
            bar.set_facecolor(CMAP(norm(temps[i])))

    #ax.bar(n, bins, histtype='bar', rwidth=0.95, stacked = True)
    title = str(Etype) + " n Values"
    ax.set_title(title)
    ax.set_xlabel("n Value")
    ax.set_ylabel("Count")


    colorbar_colors = ScalarMappable(norm)
    colorbar_colors.set_cmap(CMAP)
    fig.colorbar(colorbar_colors, ax=ax, label = c_label)


def generate_plot(ax, Eij, ex, t, slopes, use_Eij = True, cutoff = True, balanced_average = True):
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


    
    data = ax.scatter(x_, y_, label=str(ex.temp) + " °C", c=temps, s = 2, norm=colors.Normalize(MIN_TEMP, MAX_TEMP), cmap = CMAP_TEMP)

    ax.set_xlabel('Target Change')


    ax.set_xlabel("log(Tau)")

    ax.grid()
    #ax.legend(fontsize=7)
    ax.set_xscale('log')
    ax.set_yscale('log')
    if balanced_average == True:
        slopes[ax] = pd.concat([slopes[ax], get_balanced_average(np.log(x_[1:]), np.log(y_[1:]), temps)], ignore_index=True)
        return data
    elif balanced_average == "individual":
        slopes[ax] = pd.concat([slopes[ax], get_individual_slopes(np.log(x_[1:]), np.log(y_[1:]), temps)], ignore_index=True)
        return data
    else:
        m, b = np.polyfit(np.log(x_[1:]), np.log(y_[1:]), 1)
    slopes[ax].append(m)

    return data


scope = []

EXP_TYPE = EXP

for x in EXP_TYPE:
    for tem in TEMPS:
        
        if x != "ss" :
            e1 = GeneratedExperiment(x, "zz", temp = tem) 
        else:
            e1 = GeneratedExperiment(x, "xz", temp = tem) 

        print(f"GAMMA: {e1.Gamma}, LAMD: 0.15, TEMP: {e1.temp}, EXP: {e1.exptype}")
        scope.append(e1)


    plot_all_enhancement(scope, strain_over_time, balanced_average="individual", cutoff=False, fname=x)
    scope = []

