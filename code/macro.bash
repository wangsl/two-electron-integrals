#!/bin/bash

k=0
for((i=1; i<=12; i++)); do
  for((j=0; j<=12; j++)); do
    printf " %-82s // %3d \n" "i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_u_${i}_${j});" $k
    k=$((k+1))
  done
done
