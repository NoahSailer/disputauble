import numpy as np
from velocileptors.LPT.lpt_rsd_fftw import LPT_RSD
from gauss import GaussLike

# Observables: 
#      P0(k), P2(k), P4(k), Ckg
# Need to specify redshifts, fid cosmo, fid nuisance and fitting range: 
#      z, n(z), b1(z), b2(z), bs(z), SN2(z), alpha0(z), 
#      alpha2(z), alphaX(z), kmin(z), kmax(z) 








def get_closest(p):
##############################################################################################################
# My homework
    # extract LPT tables and save them in the dictionary below (for every redshift???)
    # add real space
    # make a module for adding Ckg (or add the ingredients for them)
##############################################################################################################
    # The loss function should look something like this
    # like = GaussLike(dat, cov, tmp_priors=None)
    # return like.bfchi2(thy)
    #
    lpt = LPT_RSD(k,Pk[0,:]*h**3,kIR=0.2)
    lpt.make_pltable(f[0],kmin=min(k),kmax=max(k),nk=len(k),nmax=4,apar=1,aperp=1)
    p0ktable,p2ktable,p4ktable = lpt.p0ktable/h**3,lpt.p2ktable/h**3,lpt.p4ktable/h**3