#!/bin/bash -l
#SBATCH -J add_data
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -o log/add_data.out
#SBATCH -e log/add_data.err
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH -A desi

adddata=$1
olddata=$2
conda activate cobaya_up2d8
export OMP_NUM_THREADS=8
srun -n 1  -c 8 python add_data.py ${adddata} ${olddata}
srun -n 16 -c 8 cobaya-run --force ${adddata}_${olddata}.yaml
srun -n 1  -c 8 rm ${adddata}_${olddata}.yaml