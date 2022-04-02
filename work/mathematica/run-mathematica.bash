#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mem=40GB
#SBATCH --job-name=rys
#SBATCH --array=0-155

module purge
module load mathematica/13.0.0

found=0
k=0
for((np=1; np<=12; np++)); do
  for((a=0; a<=12; a++)); do
    if [ $k -eq $SLURM_ARRAY_TASK_ID ]; then 
      found=1
      break
    else
      k=$((k+1))
    fi
  done
  if [ $found -eq 1 ]; then break; fi
done

echo 
echo "Parameters: np=$np a=$a"
FitU=1
if [ $FitU -eq 0 ]; then  
  echo "*** To fit Rys roots"
else
  echo "*** To fit Rys U"
fi
echo 

cat<< EOF | math

PrependTo[\$Path, "/home/wang/rys/two-electron-integrals/mathematica"];

<<RysChebyshevFitting\`

np = $np;
m = If[np == 1 || np == 2, 30, 25];

a = $a*10;
b = a+10;

Print["np=", np, " m=", m, " a=", a, " b=", b]

RysChebyshevFitting[np, a, b, m, $FitU];

EOF
