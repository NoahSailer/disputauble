import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager
import matplotlib.patheffects as path_effects
from getdist.mcsamples import MCSamplesFromCobaya
from getdist.mcsamples import loadMCSamples
import getdist.plots as gdplt
import yaml
from yaml.loader import SafeLoader
from scipy.stats import chi2
from scipy.special import erf
from scipy.optimize import fsolve
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
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath, amssymb}'

RM1_CUT = 0.02

################################################
####          Helper functions               ###
################################################

def load_chains(names,rm1_cut=RM1_CUT,verbose=True):
    chains = []
    for name in names:
        try: 
            chain = loadMCSamples(f'chains/{name}',settings={'ignore_rows':0.3})
            chains.append(chain)
        except: 
            s = f"Chains have not been run for {name}.yaml"
            raise RuntimeError(s)
        rm1 = chain.getGelmanRubin()
        if verbose: print(f"R-1={rm1:0.3f} for {name}.yaml",flush=True)
        if rm1>rm1_cut: 
            s = f"R-1={rm1:0.3f}>{rm1_cut} for {name}.yaml"
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

def deltaChi2_w0wa_vs_lcdm(w0wa_fn,lcdm_fn):
    w0wa_info = yaml.load(open(f'chains/{w0wa_fn}.minimize.input.yaml'), Loader=SafeLoader)
    lcdm_info = yaml.load(open(f'chains/{lcdm_fn}.minimize.input.yaml'), Loader=SafeLoader)
    w0wa_logP = -1*get_bestFit_values(w0wa_fn)['minuslogpost']
    lcdm_logP = -1*get_bestFit_values(lcdm_fn)['minuslogpost']
    w0wa_area = w0wa_info['params']['w']['prior']['max']-w0wa_info['params']['w']['prior']['min']
    w0wa_area*= w0wa_info['params']['wa']['prior']['max']-w0wa_info['params']['wa']['prior']['min']
    return 2*(w0wa_logP+np.log(w0wa_area)-lcdm_logP)

def deltaChi2_to_sigma(dchi2,dof=2):
    chi2cdf = chi2.cdf(dchi2, df=dof)
    def dcdf(x): return erf(x/np.sqrt(2)) - chi2cdf
    return fsolve(dcdf, 2)[0]

def lcdm_H0rd_OmM(fn):
    bf_H0rd = get_bestFit_values(fn)['H0rd']
    bf_OmM  = get_bestFit_values(fn)['OmegaM']
    chain   = load_chains([fn],verbose=False)[0]
    param_names_all = [p.name for p in chain.getParamNames().names]
    selected_params = ['H0rd','OmegaM']
    indices = [param_names_all.index(name) for name in selected_params]
    return np.array([bf_H0rd,bf_OmM]),chain.cov()[np.ix_(indices, indices)]

def lcdm_tension(fn1,fn2):
    bf1,cov1 = lcdm_H0rd_OmM(fn1)
    bf2,cov2 = lcdm_H0rd_OmM(fn2)
    dchi2 = np.dot(bf1-bf2,np.dot(np.linalg.inv(cov1+cov2),bf1-bf2))
    return dchi2,deltaChi2_to_sigma(dchi2,dof=2)

################################################
####               Figure 1                  ###
################################################

def make_figure_1a():
    g = gdplt.get_single_plotter(width_inch=4, ratio=8*0.63/4)
    g.settings.alpha_filled_add = 0.8
    g.settings.axes_labelsize = 28 
    g.settings.axes_fontsize = 20 
    g.settings.axis_marker_color = 'k'
    g.settings.axis_marker_lw = 1.2
    FIG0_NAMES = [
    'lcdm_mnu=0.06_tau=free_cmb-p+cmb-l',
    ]
    chains = load_chains(FIG0_NAMES)
    labels = [r'CMB']
    colors = ['green']
    g.rectangle_plot(['tau'],['OmegaM'],plot_roots=[[chains]],
                     filled=[True],colors=colors,ls=['-'],lws=[1])
    for label,c in zip(labels,colors): plt.hist([],color=c,label=label)
    ax = g.subplots[0,0]
    mu,std = 0.2967859144496603,0.00790605489981
    ax.fill_between([0.0,0.9],[mu-2*std]*2,[mu+2*std]*2,color='C0',alpha=0.3,zorder=0)
    ax.fill_between([0.0,0.9],[mu-std]*2,[mu+std]*2,color='C0',alpha=0.7,zorder=0)
    tau = np.linspace(0.01,0.2,500)
    ax.plot(tau,0.349-0.505*tau,c='yellow',ls='--',lw=4)
    plt.xlabel(r'$\tau$')
    plt.ylabel(r'$\Omega_{\rm m}$')
    legend = plt.legend(loc='upper right',frameon=False,fontsize=20)
    legend.get_frame().set_edgecolor('w')
    plt.xlim(0.03,0.12)
    plt.yticks([0.28,0.30,0.32,0.34])
    plt.ylim(0.276,0.34)
    plt.xticks([0.03,0.06,0.09,0.12])
    plt.savefig('figures/tau_OmM_contours.pdf', dpi=100, bbox_inches='tight')

def make_figure_1b():
    g = gdplt.get_single_plotter(width_inch=7, ratio=0.75)
    g.settings.alpha_filled_add = 0.8
    g.settings.axes_labelsize = 28 
    g.settings.axes_fontsize = 20 
    g.settings.axis_marker_color = 'k'
    g.settings.axis_marker_lw = 1.2
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
    ax = g.subplots[0,0]
    for i,chain in enumerate(chains):
        dist = chain.get1DDensity('H0rd'); I = np.where((dist.x>96.8)&(dist.x<103.5))
        ax.plot(dist.x[I],0.34+dist.P[I]/110,color=colors[i],ls=['-','-','--'][i],lw=3,clip_on=False)
        if i<2: ax.fill_between(dist.x[I],np.ones_like(dist.x[I])*0.34,0.34+dist.P[I]/110,
                                color=colors[i],clip_on=False,alpha=0.4)
        dist = chain.get1DDensity('OmegaM'); I = np.where((dist.x>0.276)&(dist.x<0.34))
        ax.plot(103.5+dist.P[I]/2,dist.x[I],color=colors[i],ls=['-','-','--'][i],lw=3,clip_on=False)
        if i<2: ax.fill_betweenx(dist.x[I],np.ones_like(dist.x[I])*103.5,103.5+dist.P[I]/2,
                                color=colors[i],clip_on=False,alpha=0.4)
    plt.xlabel(r'$H_0\,r_{\rm d}$ [100 km s$^{-1}$]')
    plt.ylabel('')
    plt.legend(loc='lower left',frameon=False,fontsize=22)
    plt.xticks([98,100,102])
    plt.yticks([0.28,0.30,0.32,0.34],labels=[])
    plt.xlim(96.8,103.5)
    plt.ylim(0.276,0.34)
    plt.savefig('figures/H0rd_OmM_contours.pdf', dpi=100, bbox_inches='tight')
    dchi2,nsig = lcdm_tension(FIG1_NAMES[0],FIG1_NAMES[1])
    print(f"DESI vs CMB (tau=0.06): Delta chi^2 = {dchi2:0.02f}, or {nsig:0.03f}sigma")
    dchi2,nsig = lcdm_tension(FIG1_NAMES[0],FIG1_NAMES[2])
    print(f"DESI vs CMB (tau=0.09): Delta chi^2 = {dchi2:0.02f}, or {nsig:0.03f}sigma")

def make_figure_1():
    make_figure_1a()
    make_figure_1b()

################################################
####               Figure 2                  ###
################################################

def make_figure_2():
    g = gdplt.get_single_plotter(width_inch=7, ratio=0.75)
    g.settings.alpha_filled_add = 0.4
    g.settings.axes_labelsize = 28 
    g.settings.axes_fontsize = 20 
    g.settings.axis_marker_color = 'k'
    g.settings.axis_marker_lw = 1.2
    FIG2_NAMES = [
    'lcdm_mnu=free_tau=0.06_cmb-p+cmb-l+bao', 
    'lcdm_mnu=free_tau=0.09_cmb-p+cmb-l+bao', 
    'lcdm_mnu=free_tau=free_cmb-p+cmb-lowl+cmb-l+bao'
    ]
    chains = load_chains(FIG2_NAMES)
    for c in chains:
        for x in c.getParamNames().names:
            if x.name == 'mnu_eff':
                x.label = r"M_{\nu,\mathrm{eff}}" 
    
    labels = [r'$\tau=0.06$', 
              r'$\tau=0.09$', 
              r'w/ low-$\ell$ CMB']
    colors = ['C6','k','tab:orange']
    linestyles = ['solid', 'solid', 'dashed']
    g.plot_1d(chains, 'mnu_eff', filled=True, colors=colors, ls=linestyles)
    for label,c,ls in zip(labels,colors,linestyles): plt.plot([],color=c,label=label,linestyle=ls)
    plt.legend(loc='upper left', ncol=3, fontsize=15, frameon=True, title=r'$\Lambda$CDM$ + M_{\nu,\mathrm{eff}}$', title_fontsize=20)
    plt.axvspan(-0.5, 0.0, color='grey', fill=True, alpha=0.2)
    
    plt.axvline(0.059, color='k', linestyle='dotted')
    plt.text(0.064, 0.8, "normal hierarchy\n minimum mass", fontsize=15)
    
    plt.ylim(0, 1.32)
    plt.xlim(-0.28, 0.18)
    plt.savefig("figures/mnu_eff.pdf", dpi=100, bbox_inches='tight')

################################################
####               Figure 3                  ###
################################################

def make_figure_3():
    g = gdplt.get_single_plotter(width_inch=8, ratio=0.6)
    g.settings.alpha_filled_add = 0.4
    g.settings.axes_labelsize = 28 
    g.settings.axes_fontsize = 20 
    g.settings.axis_marker_color = 'k'
    g.settings.axis_marker_lw = 1.2
    lcdm_color='gold'
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
    plt.axhline(y=0,c=lcdm_color,lw=1)
    plt.axvline(x=-1,c=lcdm_color,lw=1)
    plt.axhline(y=0,c='k',lw=2,ls='dotted')
    plt.axvline(x=-1,c='k',lw=2,ls='dotted')
    plt.xlabel(r'$w_0$')
    plt.ylabel(r'$w_a$')
    legend = plt.legend(loc=(0.6,0.55),frameon=True,fontsize=20,
                        framealpha=0.9,title=r'DESI+CMB',title_fontsize=30)
    legend.get_frame().set_edgecolor('w')
    plt.xticks([-1,-0.5,0])
    plt.yticks([-3,-2,-1,0,1])
    plt.xlim(-1.3,0.3)
    plt.ylim(-3,1)
    for i,name in enumerate(FIG3_NAMES):
        bf = get_bestFit_values(name)
        w0,wa = bf['w'],bf['wa']
        text=plt.text(w0,wa,r'$\blacktriangle$',ha='center',va='center',fontsize=15,color=['#F4C9E7','k'][i])
        text.set_path_effects([path_effects.Stroke(linewidth=2, foreground='black'),path_effects.Normal()])
    text = plt.text(-1,0,r'$\boldsymbol{\star}$',ha='center', va='center',fontsize=30,color=lcdm_color)
    text.set_path_effects([path_effects.Stroke(linewidth=2, foreground='black'),path_effects.Normal()])
    text=plt.text(-0.98,0.2,r'$\boldsymbol{\Lambda}$\textbf{CDM}',ha='left', va='bottom',fontsize=20,color=lcdm_color)
    text.set_path_effects([path_effects.Stroke(linewidth=2, foreground='black'),path_effects.Normal()])
    plt.savefig('figures/w0_wa_contours.pdf', dpi=100, bbox_inches='tight')
    for dataset in ['cmb-p+cmb-l+bao']:
        for tau in [0.06,0.09]: 
            w0wa_fn = f'w0wa_mnu=0.06_tau={tau}_{dataset}'
            lcdm_fn = f'lcdm_mnu=0.06_tau={tau}_{dataset}'
            delt = deltaChi2_w0wa_vs_lcdm(w0wa_fn,lcdm_fn)
            sig = deltaChi2_to_sigma(delt,dof=2)
            s = f"w0wa preference over lcdm -- {dataset}, tau={tau:0.03f}"
            s=s+f" -- Delta chi^2 = {delt:0.02f}, or {sig:0.03f}sigma"
            print(s,flush=True)

################################################
####               Figure 4                  ###
################################################

def make_figure_4():
    g = gdplt.get_single_plotter(width_inch=7, ratio=0.75)
    g.settings.alpha_filled_add = 0.4
    g.settings.axes_labelsize = 28 
    g.settings.axes_fontsize = 20 
    g.settings.axis_marker_color = 'k'
    g.settings.axis_marker_lw = 1.2
    FIG4_NAMES = ['lcdm_mnu=0.06_tau=free_cmb-p+cmb-l+bao',
                  'lcdm_mnu>0.06_tau=free_cmb-p+cmb-l+bao',
                  'lcdm_mnu=0.06_tau=free_cmb-p+cmb-lowl+cmb-l+bao', 
                  'lcdm_mnu>0.06_tau=free_cmb-p+cmb-lowl+cmb-l+bao'
    ]
    chains = load_chains(FIG4_NAMES)
    labels = [r"high-$\ell$ CMB + CMB lensing + BAO",
              None,
              r"+ low-$\ell$ CMB",
              None]
    colors = ['tab:purple','tab:purple', 'tab:orange', 'tab:orange']
    linestyles = ['solid', 'dashed', 'solid', 'dashed']
    g.plot_1d(chains, 'tau', colors=colors, ls=linestyles)

    for label,c,ls in zip(labels[::2],colors[::2],linestyles[::2]): plt.plot([],color=c,label=label,linestyle=ls)
    plt.plot([], color='grey', linestyle='solid', label=r"w/ $\sum m_\nu = 0.06$")
    plt.plot([], color='grey', linestyle='dashed', label=r"w/ $\sum m_\nu \geq 0.06$")
    
    plt.legend(loc='upper left', ncol=2, fontsize=14, frameon=True)
    
    plt.ylim(0, 1.27)
    plt.xlim(0.04, 0.14)
    plt.savefig("figures/tau_infer.pdf", dpi=100, bbox_inches='tight')

################################################
####          Make the figures               ###
################################################

if __name__ == "__main__":
    make_figure_1()
    make_figure_2()
    make_figure_3()
    make_figure_4()