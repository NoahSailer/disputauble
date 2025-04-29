#!/bin/bash
submit=true
njobs=3
while [[ $# -gt 0 ]]; do
  case $1 in
    --submit=*) submit="${1#*=}" ;;
    --njobs=*)  njobs="${1#*=}" ;;
    --submit)   submit="${2,,}"; shift ;;
    --njobs)    njobs="$2"; shift ;;
  esac
  shift
done
# create directories for running chains and making figures
mkdir chains log figures
# create an up-to-date conda enviornment for cobaya if one doesn't exist
if conda info --envs | awk '{print $1}' | grep -q '^cobaya_up2d8$'; then
  echo "Skipping create_cobaya_env.sh. Found existing env: cobaya_up2d8"
else
  echo "Creating conda env: cobaya_up2d8"
  bash create_cobaya_env.sh
fi
# submits jobs to the queue (if not already submitted)
cd yamls
sample_and_minimize=(
'lcdm_mnu=0.06_tau=0.06_cmb-p+cmb-l'
'lcdm_mnu=0.06_tau=0.09_cmb-p+cmb-l'
'w0wa_mnu=0.06_tau=0.06_cmb-p+cmb-l+bao'
'w0wa_mnu=0.06_tau=0.09_cmb-p+cmb-l+bao'
'lcdm_mnu=free_tau=0.06_cmb-p+cmb-l+bao'
'lcdm_mnu=free_tau=0.09_cmb-p+cmb-l+bao'
'lcdm_mnu=0.06_tau=free_cmb-p+cmb-l+bao'
)
sample=(
'lcdm_mnu=0.06_tau=free_cmb-p'
'lcdm_mnu=0.06_tau=free_cmb-p+cmb-l'
'lcdm-lite_mnu=0.06_tau=0.06_bao'
'lcdm_mnu=free_tau=free_cmb-p+cmb-lowl+cmb-l+bao'
'lcdm_mnu>0.06_tau=free_cmb-p+cmb-l+bao'
'lcdm_mnu=0.06_tau=free_cmb-p+cmb-lowl+cmb-l+bao'
'lcdm_mnu>0.06_tau=free_cmb-p+cmb-lowl+cmb-l+bao'
"${sample_and_minimize[@]}"
)
minimize=(
'lcdm_mnu=0.06_tau=0.06_cmb-p+cmb-l+bao'
'lcdm_mnu=0.06_tau=0.09_cmb-p+cmb-l+bao'
"${sample_and_minimize[@]}"
)
if $submit; then
  for filename in "${sample[@]}"; do
    if [ ! -f ../chains/$filename.chains_submitted ]; then
      JOBID=$(sbatch run_chains.sh $filename | awk '{print $4}')
      for i in $(seq 2 $njobs); do
        JOBID=$(sbatch --dependency=afterany:$JOBID run_chains.sh $filename | awk '{print $4}')
      done
      touch ../chains/$filename.chains_submitted
      echo "Submitted jobs: run_chains.sh $filename"
    else
      echo "Already submitted: run_chains.sh $filename"
    fi
  done
  for filename in "${minimize[@]}"; do
    if [ ! -f ../chains/$filename.minimizer_submitted ]; then
      JOBID=$(sbatch minimize.sh $filename | awk '{print $4}')
      touch ../chains/$filename.minimizer_submitted
      echo "Submitted job: minimize.sh $filename"
    else
      echo "Already submitted: minimize.sh $filename"
    fi
  done
fi
# make plots
cd ..
module load texlive
module load python/3.9
conda activate cobaya_up2d8
python do_analysis.py