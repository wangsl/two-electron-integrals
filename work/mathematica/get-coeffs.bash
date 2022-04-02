#!/bin/bash

n=5;

for((i=0; i<10; i++)); do
  log_file="../mathematica/slurm-15364339_$i.out"

  roots_file="rys-coeffs-$n-$i.txt"
  weights_file="rys-weights-$n-$i.txt"

  {
    grep "np=" $log_file | awk '{printf "// %-40s\n", $0}'
    egrep "^Coeffs: " $log_file | awk '{printf "%+40s, // %3d%3d%3d %4d\n", $5, $2, $3, $4, NR}'
  } > $roots_file

  {
    grep "np=" $log_file | awk '{printf "// %-40s\n", $0}'
    egrep "^Coeffs: " $log_file | awk '{printf "%+40s, // %3d%3d%3d %4d\n", $6, $2, $3, $4, NR}'
  } > $weights_file

done
