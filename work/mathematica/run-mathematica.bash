#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mem=40GB
#SBATCH --job-name=math
#SBATCH --array=10-129

module purge
module load mathematica/13.0.0

np=$((SLURM_ARRAY_TASK_ID/10))
a=$((SLURM_ARRAY_TASK_ID%10))

echo "Parameters: np=$np a=$a"

cat<< EOF | math

PrependTo[\$Path, "/home/wang/rys/two-electron-integrals/mathematica"];

<<RysChebyshevFitting\`

np = $np;
m = If[np == 1 || np == 2, 30, 25];

a = $a*10;
b = a+10;

Print["np=", np, " m=", m, " a=", a, " b=", b]

RysChebyshevFitting[np, a, b, m, 1];

EOF
