#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mem=40GB
#SBATCH --job-name=math
#SBATCH --array=0-9

module purge
module load mathematica/13.0.0

a=$SLURM_ARRAY_TASK_ID

cat<< EOF | math
<<RysChebyshevFitting\`

np = 2;
m = 30;

a = $a*10;
b = a+10;

Print["np=", np, " m=", m, " a=", a, " b=", b]

RysChebyshevFitting[np, a, b, m];
EOF
