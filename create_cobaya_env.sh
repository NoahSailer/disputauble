#!/bin/bash
# based on https://github.com/martinjameswhite/CobayaLSS/blob/main/create_cobaya_env_perlmutter.sh
module load python
conda create --name cobaya_up2d8 --clone nersc-mpi4py
conda activate cobaya_up2d8
conda install -c conda-forge numpy scipy matplotlib -y
conda install -c conda-forge astropy sympy pandas cython -y
conda install -c conda-forge pyfftw healpy -y
conda install -c conda-forge ipykernel ipython jupyter -y
python3 -m ipykernel install --user --name cobaya_up2d8 --display-name cobaya-up2d8-env
python3 -m pip install cobaya --upgrade
if [ ! -d $SCRATCH/Cobaya ]; then
  mkdir $SCRATCH/Cobaya
else
  rm -rf $SCRATCH/Cobaya/Packages
fi
cobaya-install cosmo -p $SCRATCH/Cobaya/Packages
cobaya-install cosmo --upgrade -p $SCRATCH/Cobaya/Packages
python3 -m pip install act_dr6_lenslike  --user
cobaya-install act_dr6_lenslike.ACTDR6LensLike -p $SCRATCH/Cobaya/Packages
python3 -m pip install -v git+https://github.com/gerrfarr/CMBLensLklh --user
cobaya-install CMBLensLklh.CMBLensLklh -p $SCRATCH/Cobaya/Packages