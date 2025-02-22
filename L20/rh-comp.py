from lagrange_gen import Experiment
import matplotlib.pyplot as plt
import numpy as np
import os

import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredText



T = 10e-5


def plot_experiment_enhancement(ex, T):
    '''
    ex is an Experiment object which already has generated an .nc file

    '''
    df = ex.get_dataframe()



    dpi, scale = 200, 3.3
    fig = plt.figure(figsize=(3/2*2.1*scale,2.3*scale/1.7))
    gs = gridspec.GridSpec(1,2, height_ratios=[.5], width_ratios=[1,1])
    #gs.update(left=-0.03, right=1-0.06/3, top=0.97, bottom=0.20, wspace=0.015*18, hspace=0.35)

    ax_Elin    = plt.subplot(gs[0, 0])
    ax_Enlin   = plt.subplot(gs[0, 1])
  
    lw0, lw1, lw2 = 2.5,2.25,2.0


    lblm = lambda ii,jj: '$E_{m_%i m_%i}$'%(ii+1,jj+1)
    lblp = lambda ii,jj: '$E_{p_%i p_%i}$'%(ii+1,jj+1)

    def plot_enhancements(ax, Eij):
        df['tau'] = Eij / T
        changes = [df.tau[1] - df.tau[0]] * len(df.step)

        
        ax.semilogy(changes, Eij, '-', c="red",   label=lblm(0,0), lw=lw0)

        #xlims=[0, df.step[len(df.step) - 1]]
        #ax.semilogy(xlims, [2.5,2.5], '--k', lw=1) 
        #ax.semilogy(xlims, [4.375,4.375], '--k', lw=1)
        #ax.set_xlim(xlims)
        #ax.set_ylim(np.amin(Eij[:]), np.amax(Eij[:]))

        
        if ex.exptype == "ss":
            ax.set_xlabel('Target Angle')
        else:
            ax.set_xlabel('Target Change')

        ax.set_xlabel("ė")

        ax.set_ylabel('$E_{vw}$')
        ax.grid()

    plot_enhancements(ax_Enlin, df.nonlinear_enhancement)
    ax_Enlin.set_title(r'Nonlinear Sachs')

    plot_enhancements(ax_Elin,  df.linear_enhancement)
    ax_Elin.set_title(r'Linear mixed Taylor--Sachs')



    os.makedirs("output", exist_ok=True)
    fout = 'output/cache.png'
    print('Saving %s'%(fout))
    plt.savefig(fout, dpi=dpi)
    plt.close('all')



e1 = Experiment("ss", "xz", temp = -30) 

plot_experiment_enhancement(e1, T)