#!/bin/bash
reanalyze=true
continue_chains=false
njobs=4
while [[ "$#" -gt 0 ]]; do
  [[ "$1" == "--reanalyze" ]] && reanalyze="${2,,}" && shift
  shift
done
# create directories for running chains and making figures
mkdir chains log figures
# create an up to date conda enviornment for cobaya if one doesn't exist
if conda info --envs | awk '{print $1}' | grep -q '^cobaya_up2d8$'; then
  echo "Skipping create_cobaya_env.sh. Found existing env: cobaya_up2d8"
else
  echo "Creating conda env: cobaya_up2d8"
  bash create_cobaya_env.sh
fi
# !!CURRENTLY RESTRICTED TO FIGS 1 AND 3!!
cd yamls
if $reanalyze; then
  for filename in *tau\=0.0*.yaml; do
  # If not already submitted:
  #   - Submits njobs run_chains jobs for each yaml file. 
  #     The njobs are daisy-chained by their dependencies.
  #   - Submits minimizer job for each yaml file.
    base="${filename%.yaml}"
    # chains
    if [ ! -f ../chains/$base.chains_submitted ]; then
      JOBID=$(sbatch run_chains.sh $base | awk '{print $4}')
      for i in $(seq 2 $njobs); do
        JOBID=$(sbatch --dependency=afterany:$JOBID run_chains.sh $base | awk '{print $4}')
      done
      touch ../chains/$base.chains_submitted
      echo "Submitted jobs: run_chains.sh $base"
    else
      echo "Already submitted: run_chains.sh $base"
    fi
    # minimizer
    if [ ! -f ../chains/$base.minimizer_submitted ]; then
      JOBID=$(sbatch minimize.sh $base | awk '{print $4}')
      touch ../chains/$base.minimizer_submitted
      echo "Submitted job: minimize.sh $base"
    else
      echo "Already submitted: minimize.sh $base"
    fi
  done
fi
# make plots
cd ..
module load texlive
module load python/3.9
conda activate cobaya_up2d8
python do_analysis.py