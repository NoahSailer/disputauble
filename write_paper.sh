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
# run chains & minimizer (CURRENTLY RESTRICTED TO FIG 1 AND 3)
cd yamls
for filename in *tau\=0.0*.yaml; do
  base="${filename%.yaml}"
  if [ ! -f ../chains/$base.chains_submitted ]; then
    echo "Submitted job: run_chains.sh $base"
    sbatch run_chains.sh $base
    touch ../chains/$base.chains_submitted
  else
    echo "Already submitted: run_chains.sh $base"
  fi
  if [ ! -f ../chains/$base.minimizer_submitted ]; then
    echo "Submitted job: minimize.sh $base"
    sbatch minimize.sh $base
    touch ../chains/$base.minimizer_submitted
  else
    echo "Already submitted: minimize.sh $base"
  fi
done
# make plots
cd ..
module load texlive
module load python/3.9
conda activate cobaya_up2d8
python make_plots.py
