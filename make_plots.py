import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager
from getdist.mcsamples import MCSamplesFromCobaya
from getdist.mcsamples import loadMCSamples
import getdist.plots as gdplt
# make plots pretty
plt.rc('font',**{'size':'22','family':'serif','serif':['CMU serif']})
plt.rc('mathtext', **{'fontset':'cm'})
plt.rc('text', usetex=True)
plt.rc('legend',**{'fontsize':'18'})
matplotlib.use(plt.rcParams["backend"])
matplotlib.rcParams['axes.linewidth'] = 3
matplotlib.rcParams['axes.labelsize'] = 30
matplotlib.rcParams['xtick.labelsize'] = 25 
matplotlib.rcParams['ytick.labelsize'] = 25
matplotlib.rcParams['legend.fontsize'] = 25
matplotlib.rcParams['xtick.major.size'] = 10
matplotlib.rcParams['ytick.major.size'] = 10
matplotlib.rcParams['xtick.minor.size'] = 5
matplotlib.rcParams['ytick.minor.size'] = 5
matplotlib.rcParams['xtick.major.width'] = 3
matplotlib.rcParams['ytick.major.width'] = 3
matplotlib.rcParams['xtick.minor.width'] = 1.5
matplotlib.rcParams['ytick.minor.width'] = 1.5
matplotlib.rcParams['axes.titlesize'] = 30
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

RM1_CUT = 1e5

################################################
####   Check if chains have been analyzed    ###
################################################

def load_chains(names,rm1_cut=RM1_CUT):
    chains = []
    for name in names:
        try: 
            chain = loadMCSamples(f'chains/{name}',settings={'ignore_rows':0.3})
            chains.append(chain)
        except: 
            s = f"make_plots.py failed. Chains have not been run for {name}.yaml"
            raise RuntimeError(s)
        rm1 = chain.getGelmanRubin()
        if rm1>rm1_cut: 
            s = f"make_plots.py failed. R-1={rm1:0.3f}>{rm1_cut} for {name}.yaml"
            raise RuntimeError(s)
    return chains

def load_minima(names):
    minima = []
    for name in names: 
        continue
    return minima

################################################
####               Figure 1                  ###
################################################

def make_figure_1():
    g = gdplt.get_single_plotter(width_inch=7, ratio=0.8)
    g.settings.alpha_filled_add = 0.8
    g.settings.axes_labelsize = 28 
    g.settings.axes_fontsize = 20 
    g.settings.axis_marker_color = 'k'
    g.settings.axis_marker_lw = 1.2
    g.settings.figure_legend_ncol = 2
    g.settings.linewidth_contour = 2
    FIG1_NAMES = [
    'lcdm-lite_mnu=0.06_tau=0.06_bao',
    'lcdm_mnu=0.06_tau=0.06_cmb-p+cmb-l',
    'lcdm_mnu=0.06_tau=0.09_cmb-p+cmb-l',
    ]
    chains = load_chains(FIG1_NAMES)
    labels = [r'DESI DR2 BAO',r'CMB, $\tau=0.06$',r'CMB, $\tau=0.09$']
    colors = ['C0','C6','k']
    g.rectangle_plot(['H0rd'],['OmegaM'],plot_roots=[[chains]],
                     filled=[True,True,False],colors=colors,ls=['-','-','--'])
    for label,c in zip(labels,colors): plt.hist([],color=c,label=label)
    plt.xlabel(r'$H_0\,r_{\rm d}$ [100 km s$^{-1}$]')
    plt.ylabel(r'$\Omega_{\rm m}$')
    plt.legend(loc='lower left',frameon=False,fontsize=20)
    plt.xticks([98,100,102])
    plt.yticks([0.28,0.30,0.32,0.34])
    plt.xlim(96.8,103.5)
    plt.ylim(0.276,0.34)
    plt.savefig('figures/H0rd_OmM_contours.pdf', dpi=100, bbox_inches='tight')

################################################
####               Figure 3                  ###
################################################

def make_figure_3():
    g = gdplt.get_single_plotter(width_inch=7, ratio=0.8)
    g.settings.alpha_filled_add = 0.8
    g.settings.axes_labelsize = 28 
    g.settings.axes_fontsize = 20 
    g.settings.axis_marker_color = 'k'
    g.settings.axis_marker_lw = 1.2
    g.settings.figure_legend_ncol = 2
    g.settings.linewidth_contour = 2
    FIG3_NAMES = [
    'w0wa_mnu=0.06_tau=0.06_cmb-p+cmb-l+bao',
    'w0wa_mnu=0.06_tau=0.09_cmb-p+cmb-l+bao',
    ]
    chains = load_chains(FIG3_NAMES)
    labels = [r'$\tau=0.06$',r'$\tau=0.09$']
    colors = ['C0','C3']
    g.rectangle_plot(['w'],['wa'],plot_roots=[[chains]],
                     filled=[True,True],colors=colors,ls=['-','-'])
    for label,c in zip(labels,colors): plt.hist([],color=c,label=label)
    plt.text(-0.32,-0.5,r'DESI+CMB',fontsize=30)
    plt.axhline(y=0,c='k',lw=2,ls='--')
    plt.axvline(x=-1,c='k',lw=2,ls='--')
    plt.xlabel(r'$w_0$')
    plt.ylabel(r'$w_a$')
    plt.legend(loc=(0.6,0.55),frameon=False,fontsize=20)
    plt.xticks([-1,-0.5,0])
    plt.yticks([-3,-2,-1,0])
    plt.xlim(-1.1,0.3)
    plt.ylim(-3,0.3)
    plt.text(-1,0,r'$\boldsymbol{\star}$',ha='center', va='center',fontsize=60,color='goldenrod')
    plt.text(-0.98,-0.15,r'$\boldsymbol{\Lambda}$\textbf{CDM}',ha='left', va='top',fontsize=30,color='goldenrod')
    plt.savefig('figures/w0_wa_contours.pdf', dpi=100, bbox_inches='tight')

################################################
####          Make the figures               ###
################################################

if __name__ == "__main__":
    make_figure_1()
    make_figure_3()