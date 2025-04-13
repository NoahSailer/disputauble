import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager
from getdist.mcsamples import MCSamplesFromCobaya
from getdist.mcsamples import loadMCSamples
import getdist.plots as gdplt
import yaml
from yaml.loader import SafeLoader
from scipy.stats import chi2, norm
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
####          Helper functions               ###
################################################

def load_chains(names,rm1_cut=RM1_CUT,verbose=True):
    chains = []
    for name in names:
        try: 
            print(name)
            chain = loadMCSamples(f'chains/{name}',settings={'ignore_rows':0.3})
            chains.append(chain)
        except: 
            s = f"make_plots.py failed. Chains have not been run for {name}.yaml"
            raise RuntimeError(s)
        rm1 = chain.getGelmanRubin()
        if verbose: print(f"R-1={rm1:0.3f} for {name}.yaml",flush=True)
        if rm1>rm1_cut: 
            s = f"make_plots.py failed. R-1={rm1:0.3f}>{rm1_cut} for {name}.yaml"
            raise RuntimeError(s)
    return chains

def get_bestFit_values(fn):
    fullfn = f'chains/{fn}.minimum.txt'
    for line in open(fullfn):
        li=line.strip()
        if li.startswith("#"):
            header = li
    header = np.array(header.split())[1:].tolist()
    bf = np.loadtxt(fullfn)
    return dict(zip(header,bf))

################################################
####            Tension metrics              ###
################################################

def deltaChi2_w0wa_vs_lcdm(w0wa_fn,lcdm_fn):
    w0wa_info = yaml.load(open(f'chains/{w0wa_fn}.minimize.input.yaml'), Loader=SafeLoader)
    lcdm_info = yaml.load(open(f'chains/{lcdm_fn}.minimize.input.yaml'), Loader=SafeLoader)
    w0wa_logP = -1*get_bestFit_values(w0wa_fn)['minuslogpost']
    lcdm_logP = -1*get_bestFit_values(lcdm_fn)['minuslogpost']
    w0wa_area = w0wa_info['params']['w']['prior']['max']-w0wa_info['params']['w']['prior']['min']
    w0wa_area*= w0wa_info['params']['wa']['prior']['max']-w0wa_info['params']['wa']['prior']['min']
    return 2*(w0wa_logP+np.log(w0wa_area)-lcdm_logP)

def deltaChi2_w0wa_vs_lcdm_tau(tau,dataset='cmb-p+cmb-l+bao'):
    w0wa_fn = f'w0wa_mnu=0.06_tau={tau:0.02f}_{dataset}'
    lcdm_fn = f'lcdm_mnu=0.06_tau={tau:0.02f}_{dataset}'
    return deltaChi2_w0wa_vs_lcdm(w0wa_fn,lcdm_fn)

def deltaChi2_to_sigma(dchi2,dof=2):
    chi2cdf = chi2.cdf(dchi2, df=dof)
    x = np.linspace(-10,10,50000)
    normcdf = norm.cdf(x)
    return x[np.argmin((chi2cdf-normcdf)**2)]

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
    colors = ['C6','k']
    g.rectangle_plot(['w'],['wa'],plot_roots=[[chains]],
                     filled=[True,False],colors=colors,ls=['-','--'])
    for label,c in zip(labels,colors): plt.hist([],color=c,label=label)
    plt.text(-0.5,0.2,r'DESI+CMB',fontsize=30)
    plt.axhline(y=0,c='gray',lw=2,ls='dotted')
    plt.axvline(x=-1,c='gray',lw=2,ls='dotted')
    plt.xlabel(r'$w_0$')
    plt.ylabel(r'$w_a$')
    plt.legend(loc=(0.6,0.55),frameon=False,fontsize=20)
    plt.xticks([-1.5,-1,-0.5,0])
    plt.yticks([-3,-2,-1,0,1])
    plt.xlim(-1.5,0.3)
    plt.ylim(-3,1)
    for i,name in enumerate(FIG3_NAMES):
        bf = get_bestFit_values(name)
        w0,wa = bf['w'],bf['wa']
        plt.text(w0,wa,r'$\boldsymbol{\spadesuit}$',ha='center', va='center',fontsize=15,color=colors[i])
    plt.text(-1,0,r'$\boldsymbol{\star}$',ha='center', va='center',fontsize=30,color='gray')
    plt.text(-0.98,0.2,r'$\boldsymbol{\Lambda}$\textbf{CDM}',ha='left', va='bottom',fontsize=20,color='gray')
    plt.savefig('figures/w0_wa_contours.pdf', dpi=100, bbox_inches='tight')

################################################
####          Make the figures               ###
####        Print tension metrics            ###
################################################

if __name__ == "__main__":
    make_figure_1()
    make_figure_3()
    for dataset in ['cmb-p+cmb-l+bao']:
        for tau in [0.06,0.09]: 
            delt = deltaChi2_w0wa_vs_lcdm_tau(tau=tau,dataset=dataset)
            sig = deltaChi2_to_sigma(delt,dof=2)
            s = f"w0wa preference over lcdm -- {dataset}, tau={tau:0.02f}"
            s=s+f" -- Delta chi^2 = {delt:0.02f}, or {sig:0.01f}sigma"
            print(s,flush=True)