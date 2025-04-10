#!/bin/bash -l
#SBATCH -J minimize
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -o log/minimize.out
#SBATCH -e log/minimize.err
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH -A desi

name=$1
conda activate cobaya_up2d8
rm chains/*lock*
export COBAYA_USE_FILE_LOCKING=false
export OMP_NUM_THREADS=4
srun -N 1 -n 32 -c 4 cobaya-run --minimize --force ${name}.yaml