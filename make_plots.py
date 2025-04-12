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

################################################
####               Figure 1                  ###
################################################

g = gdplt.get_single_plotter(width_inch=7, ratio=0.8)
g.settings.alpha_filled_add = 0.8
g.settings.axes_labelsize = 28 
g.settings.axes_fontsize = 20 
g.settings.axis_marker_color = 'k'
g.settings.axis_marker_lw = 1.2
g.settings.figure_legend_ncol = 2
g.settings.linewidth_contour = 2

names  = ['lcdm-lite_mnu=0.06_tau=0.06_bao',
          'lcdm_mnu=0.06_tau=0.06_cmb-p+cmb-l']
chains = [loadMCSamples(f'chains/{name}',settings={'ignore_rows':0.3}) for name in names]
for chain in chains:
    chain.addDerived(chain.samples[:,chain.index['H0rd']]/100, name='H0rd/100')
labels = [r'DESI DR2 BAO',r'CMB, $\tau=0.06$']
colors = ['C0','C6']
g.rectangle_plot(['H0rd/100'],['OmegaM'],plot_roots=[[chains]],filled=True,colors=colors)
for label,c in zip(labels,colors):
    plt.hist([],color=c,label=label)
plt.xlabel(r'$H_0\,r_d$ [100 km s$^{-1}$]')
plt.ylabel(r'$\Omega_{\rm m}$')
plt.legend(loc='lower left',frameon=False,fontsize=20)
plt.xticks([98,100,102,104])
plt.yticks([0.28,0.30,0.32,0.34])
plt.savefig('figures/H0rd_OmM_contours.pdf', dpi=100, bbox_inches='tight')