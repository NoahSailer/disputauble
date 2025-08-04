# example usage: 
#  salloc --nodes 4 --qos interactive --time 02:00:00 --constraint cpu
#  conda activate cobaya_up2d8
#  srun -N 4 -n 256 -c 2 python ppd.py
from mpi4py import MPI
import os, json
import numpy as np
from getdist import MCSamples, loadMCSamples
from cobaya.yaml import yaml_load_file
from classy import Class
from scipy.optimize import minimize
from velocileptors.LPT.lpt_rsd_fftw import LPT_RSD
from gauss import GaussLike

p_fid = {'As': 2.150860e-9,'ns': 0.9625356,
        'H0': 67.02393,'mnu': 0.06,
        'tau': 0.055,'ombh2': 0.02219218,'omch2': 0.1203058,
        'b1':3,'b2':0,'bs':0,'alpha0':0,'alpha2':0,'sn':0,'sn2':0,'z_lpt':2}

def p0p2p4_to_pkmu(p0,p2,p4,mu):
    Nmu = len(mu)
    return np.repeat(p0,Nmu)+0.5*(3*mu**2-1)*np.repeat(p2,Nmu)+0.125*(35*mu**4-30*mu**2+3)*np.repeat(p4,Nmu)


def get_closest_lpt_for_snapshot_geometry(p,p_fid,kmin,kmax,Nk,Nmu,V=1e9):
    pred    = get_prediction(p)
    pred_fid= get_prediction(p_fid)
    klin    = pred['k [h/Mpc]']
    # Set up the k-mu grid 
    k = np.logspace(np.log10(kmin),np.log10(kmax),Nk)
    dk = list(k[1:]-k[:-1])
    dk.insert(0,dk[0])
    dk = np.array(dk)
    k = np.repeat(k,Nmu)
    dk = np.repeat(dk,Nmu)
    mu = np.linspace(0.,1.,Nmu)
    dmu = list(mu[1:]-mu[:-1])
    dmu.append(dmu[-1])
    dmu = np.array(dmu)
    mu = np.tile(mu,Nk)
    dmu = np.tile(dmu,Nk)
    # make tables for sampled cosmology
    zidx    = np.argmin((pred['z']-p_fid['z_lpt'])**2)
    # FIX THIS 
    apar    = 1
    aperp   = 1
    Pcb     = LPT_RSD(klin,pred['Pcb_lin [(Mpc/h)^3]'][zidx,:],kIR=0.2)
    print('using velocileptors')
    Pcb.make_pltable(pred['f'][zidx],kmin=kmin,kmax=kmax,nk=Nk,nmax=4,apar=apar,aperp=aperp)
    print('done')
    # and fiducial prediction
    Pcb_fid = LPT_RSD(klin,pred_fid['Pcb_lin [(Mpc/h)^3]'][zidx,:],kIR=0.2)
    Pcb_fid.make_pltable(pred_fid['f'][zidx],kmin=kmin,kmax=kmax,nk=Nk,nmax=4,apar=1,aperp=1)
    bvec_fid = [p_fid['b1'],p_fid['b2'],p_fid['bs'],0.,p_fid['alpha0'],p_fid['alpha2'],0,0,p_fid['sn'],p_fid['sn2'],0]
    _,p0_fid,p2_fid,p4_fid = Pcb_fid.combine_bias_terms_pkell(bvec_fid)
    pkmu_fid  = p0p2p4_to_pkmu(p0_fid,p2_fid,p4_fid,mu)
    cov       = np.diag(4*np.pi**2*pkmu_fid**2/k**2/dk/dmu/V)
    tmp_priors= [[0,50],[0,50],[p_fid['sn'],p_fid['sn']*0.3],[p_fid['sn2'],p_fid['sn2']]]
    like      = GaussLike(pkmu_fid, cov, tmp_priors=tmp_priors)
    def chi2(bvec,get_prediciton=False):
        b1,b2,bs  = bvec
        pkmu_null = Pcb.combine_bias_terms_pkell([0,0,0,0,0,0,0,0,0,0,0])
        pkmu_a0   = Pcb.combine_bias_terms_pkell([0,0,0,0,1,0,0,0,0,0,0])
        pkmu_a2   = Pcb.combine_bias_terms_pkell([0,0,0,0,0,1,0,0,0,0,0])
        pkmu_sn   = Pcb.combine_bias_terms_pkell([0,0,0,0,0,0,0,0,1,0,0])
        pkmu_sn2  = Pcb.combine_bias_terms_pkell([0,0,0,0,0,0,0,0,0,1,0])
        pkmu_     = Pcb.combine_bias_terms_pkell([b1,b2,bs,0,0,0,0,0,0,0,0])
        pkmu = np.zeros((len(k_fid),5))
        pkmu[:,0] = p0p2p4_to_pkmu(pkmu_[1],pkmu_[2],pkmu_[3],mu) 
        pkmu[:,1] = p0p2p4_to_pkmu(pkmu_a0[1]-pkmu_null[1],pkmu_a0[2]-pkmu_null[2],pkmu_a0[3]-pkmu_null[3],mu)
        pkmu[:,2] = p0p2p4_to_pkmu(pkmu_a2[1]-pkmu_null[1],pkmu_a2[2]-pkmu_null[2],pkmu_a2[3]-pkmu_null[3],mu)
        pkmu[:,3] = p0p2p4_to_pkmu(pkmu_sn[1]-pkmu_null[1],pkmu_sn[2]-pkmu_null[2],pkmu_sn[3]-pkmu_null[3],mu)
        pkmu[:,4] = p0p2p4_to_pkmu(pkmu_sn2[1]-pkmu_null[1],pkmu_sn2[2]-pkmu_null[2],pkmu_sn2[3]-pkmu_null[3],mu)
        if get_prediciton: return like.getBestFit(pkmu)
        return like.bfchi2(pkmu)
    result = minimize(chi2, np.array([p_fid['b1'],p_fid['b2'],p_fid['bs']]))
    return {'k_lpt':k,'mu_lpt':mu,'dPkmu_lpt':chi2(result.x,get_prediciton=True)-pkmu_fid}

def convert_numpy_to_list(obj):
    """
    Recursively convert all numpy arrays in a dictionary (or list) to Python lists.
    """
    if isinstance(obj, dict):
        return {k: convert_numpy_to_list(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy_to_list(v) for v in obj]
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    else:
        return obj

def get_fixed_parameters_from_yaml(yaml_path):
    """
    Returns a dictionary of fixed parameters and their values from a Cobaya YAML file.
    """
    info = yaml_load_file(yaml_path)
    params = info.get("params", {})
    fixed_params = {}
    for name, spec in params.items():
        if isinstance(spec, (int, float)):
            fixed_params[name] = spec
        elif isinstance(spec, dict):
            if "prior" not in spec and not spec.get("derived", False):
                if "value" in spec: fixed_params[name] = spec["value"]
                else: fixed_params[name] = None
    return fixed_params

def get_sampled_parameters_from_yaml(yaml_path):
    """
    Load Cobaya YAML file and return only sampled (not fixed nor derived) parameters.
    """
    info = yaml_load_file(yaml_path)
    params = info.get("params", {})
    sampled_params = []
    for param, cfg in params.items():
        if isinstance(cfg, dict):
            if not cfg.get("derived", False) and "prior" in cfg:
                sampled_params.append(param)
        elif cfg is None:
            sampled_params.append(param)  
    return sampled_params

def get_derived_parameters_from_yaml(yaml_path):
    """
    Returns a list of parameter names that are marked as derived in a Cobaya YAML file.
    """
    info = yaml_load_file(yaml_path)
    params = info.get("params", {})
    derived_params = []
    for name, spec in params.items():
        if isinstance(spec, dict) and spec.get("derived", False):
            derived_params.append(name)
    return derived_params

def get_best_fit_values(fullfn):
    """
    Get the best fit values of all parameters.
    """
    for line in open(fullfn):
        li=line.strip()
        if li.startswith("#"):
            header = li
    header = np.array(header.split())[1:].tolist()
    bf = np.loadtxt(fullfn)
    return dict(zip(header,bf))

def load_fair_sample(yaml_root, sample_size=1000):
    """
    Load a fair sample of sampled and derived parameters from a GetDist chain.
    """
    samples = loadMCSamples(f"{yaml_root}",settings={'ignore_rows':0.3})
    params = [p.name for p in samples.getParamNames().names]
    chain_array = samples.samples
    weights = samples.weights
    weights /= np.sum(weights)
    indices = np.random.choice(len(weights), size=sample_size, p=weights)
    fair_samples = chain_array[indices]
    param_indices = [samples.index[name] for name in params]
    output = {name: fair_samples[:, i] for name, i in zip(params, param_indices)}
    return output

def get_prediction(p,can_predict=False,include_lpt=False):
    """
    Do some CLASS-y and Velocileptor-y things given input parameters (p).
    """
    output={'z':None,'k [h/Mpc]':None,'Pcb_lin [(Mpc/h)^3]':None,'Pcb_nonlin [(Mpc/h)^3]':None,
            'Pcb_nonlin_norm [(Mpc/h)^3]':None,'DH [Mpc]':None,'DM [Mpc]':None,
            'sigma8z':None,'DM_div_DH':None,'DV_div_rd':None}
    if include_lpt: output = {**output,**{'k_lpt':None,'mu_lpt':None,'dPkmu_lpt':None}}
    if can_predict: return list(output.keys())
    h = p['H0']/100
    params = {'output': 'mPk lCl','P_k_max_h/Mpc': 5.,'non linear':'hmcode', 
              'z_pk': '0.0,5.0','A_s': p['As'],'n_s': p['ns'],
              'h': h, 'N_ur': 2.0308,'N_ncdm': 1,'m_ncdm': p['mnu'],
              'tau_reio': p['tau'],'omega_b': p['ombh2'],'omega_cdm': p['omch2']}
    if 'omk' in p.keys(): 
        params['Omega_k'] = p['omk']   
    if 'wa' in p.keys():
        params['Omega_Lambda']            = 0
        params['fluid_equation_of_state'] = 'CLP'
        params['w0_fld'] = p['w']
        params['wa_fld'] = p['wa']
    cosmo = Class()
    cosmo.set(params)
    cosmo.compute()
    k = np.logspace(-3,0,250) # h/Mpc units
    z = np.linspace(1e-2,5,75)
    Hz= cosmo.Hubble(z)
    Cz= cosmo.comoving_distance(z) # Mpc units
    Dz= cosmo.scale_independent_growth_factor(z)
    s8= cosmo.sigma8()
    rd= cosmo.rs_drag()
    f = cosmo.scale_independent_growth_factor_f(z)
    Pl= np.array([[cosmo.pk_cb_lin(kk*h,zz) for kk in k] for zz in z])*h**3 # (Mpc/h)^3 units
    Pn= np.array([[cosmo.pk_cb(kk*h,zz) for kk in k] for zz in z])*h**3     # (Mpc/h)^3 units
    Pn_norm= np.array([[cosmo.pk_cb(kk*h,zz)/(s8*cosmo.scale_independent_growth_factor(zz))**2 for kk in k] for zz in z])*h**3 
    output = {'z':z,
              'k [h/Mpc]':k,
              'Pcb_lin [(Mpc/h)^3]':Pl,
              'Pcb_nonlin [(Mpc/h)^3]':Pn,
              'Pcb_nonlin_norm [(Mpc/h)^3]':Pn_norm,
              'DH [Mpc]':1/Hz,
              'DM [Mpc]':Cz,
              'f':f,
              'sigma8z':s8*Dz,
              'DM_div_DH':Cz*Hz,
              'DV_div_rd':(z*Cz**2/Hz)**(1/3)/rd
             }
    if include_lpt: output = {**output,**get_closest_lpt_for_snapshot_geometry(p,p_fid,1e-3,0.4,100,50,V=1e9)}
    return output

def ppd(yaml_root, sample_size):
    """
    Computes the PPD for quantities provided by get_prediction, in parallel using MPI.
    """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    if rank == 0:
        fixed_params = get_fixed_parameters_from_yaml(f"{yaml_root}.input.yaml")
        fair_samples = load_fair_sample(yaml_root, sample_size)
        can_predict  = get_prediction({},can_predict=True)
        try: 
            best_fit            = get_best_fit_values(f"{yaml_root}.minimum.txt")
            best_fit_prediction = get_prediction({**fixed_params, **best_fit})
        except: 
            best_fit = None
            best_fit_prediction = None
    else:
        fixed_params = None
        best_fit = None
        fair_samples = None
        best_fit_prediction = None
        can_predict = None
    fixed_params        = comm.bcast(fixed_params, root=0)
    best_fit_prediction = comm.bcast(best_fit_prediction, root=0)
    can_predict         = comm.bcast(can_predict, root=0)
    fair_samples        = comm.bcast(fair_samples, root=0)
    local_predictions   = {name: [] for name in can_predict}
    for i in range(rank, sample_size, size):
        fair_sample = {name: fair_samples[name][i] for name in fair_samples}
        prediction = get_prediction({**fixed_params, **fair_sample})
        for name in can_predict:
            local_predictions[name].append(list(prediction[name]))
    gathered_predictions = comm.gather(local_predictions, root=0)
    if rank == 0:
        fair_sample_predictions = {name: [] for name in can_predict}
        for partial in gathered_predictions:
            for name in can_predict:
                fair_sample_predictions[name].extend(partial[name])
        output = {
            'input_yaml': yaml_load_file(f"{yaml_root}.input.yaml"),
            'best_fit_prediction': best_fit_prediction,
            'mean_prediction': {name:np.mean(fair_sample_predictions[name], axis=0) for name in can_predict},
            'median_prediction': {name:np.median(fair_sample_predictions[name], axis=0) for name in can_predict},
            '2sigma_low':  {name: np.percentile(fair_sample_predictions[name], 2.5, axis=0) for name in can_predict},
            '1sigma_low':  {name: np.percentile(fair_sample_predictions[name], 16, axis=0) for name in can_predict},
            '1sigma_high': {name: np.percentile(fair_sample_predictions[name], 84, axis=0) for name in can_predict},
            '2sigma_high': {name: np.percentile(fair_sample_predictions[name], 97.5, axis=0) for name in can_predict},
            'sample_size': sample_size,
        }
        with open(f"ppd_{yaml_root.split('/')[-1]}.json", "w") as f:
            json.dump(convert_numpy_to_list(output), f, indent=2)

#if __name__ == "__main__":
#    yaml_roots   = ['../chains/lcdm_mnu=0.06_tau=0.06_cmb-p+cmb-l',
#                    '../chains/lcdm_mnu=0.06_tau=0.09_cmb-p+cmb-l+bao',
#                    '../chains/w0wa_mnu=0.06_tau=0.06_cmb-p+cmb-l+bao',
#                    '../chains/omk_mnu=0.06_tau=0.06_cmb-p+cmb-l+bao']
#    sample_size  = 1024
#    for yaml_root in yaml_roots: ppd(yaml_root, sample_size)

get_closest_lpt_for_snapshot_geometry(p_fid,p_fid,1e-3,0.4,20,10,V=1e9)