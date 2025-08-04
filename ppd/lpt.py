import numpy as np
from velocileptors.LPT.lpt_rsd_fftw import LPT_RSD
from gauss import GaussLike
from ppd import get_prediction


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

# Observables: 
#      P0(k), P2(k), P4(k), Ckg
# Need to specify redshifts, fid cosmo, fid nuisance and fitting range: 
#      z, n(z), b1(z), b2(z), bs(z), SN2(z), alpha0(z), 
#      alpha2(z), alphaX(z), kmin(z), kmax(z) 

