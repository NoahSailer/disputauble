#!/bin/bash
# create directories for running chains and making figures
mkdir chains log figures
# create an up to date conda enviornment for cobaya if one doesn't exist
if conda info --envs | awk '{print $1}' | grep -q '^cobaya_up2d8$'; then
  echo "Skipping create_cobaya_env.sh. Found existing env: cobaya_up2d8"
else
  echo "Creating conda env: cobaya_up2d8"
  bash create_cobaya_env.sh
fi
# run chains & minimizer (UPDATE)
cd yamls
for filename in *tau\=0.0*.yaml; do
  base="${filename%.yaml}"
  if [ ! -f ../chains/$base.1.txt ]; then
    echo "Submitted job: run_chains.sh $base"
    sbatch run_chains.sh $base
  else
    echo "Skipping run_chains.sh $base, found existing chains."
  fi
  if [ ! -f ../chains/$base.minimize.minumum ]; then
    echo "Submitted job: minimize.sh $base"
    sbatch minimize.sh $base
  else
    echo "Skipping minimize.sh $base, found existing minimum."
  fi
done
# make plots (ADDRESS ISSUES WITH MATPLOTLIB AT NERSC)
#cd ..
#module load texlive
#source activate noah_base
#python make_plots.py