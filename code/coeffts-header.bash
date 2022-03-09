#!/bin/bash

#job_id="15369279" # n=1
job_id="15369314" # n=2
job_id="15369155" # n=3
job_id="15369174" # n=4
job_id="15364339" # n=5
job_id="15365120" # n=6
job_id="15366995" # n=7
job_id="15367051" # n=8
job_id="15367133" # n=9
job_id="15367134" # n=10
job_id="15367157" # n=11
job_id="15367370" # n=12

for((i=0; i<10; i++)); do
  log_file="../mathematica/slurm-${job_id}_${i}.out"
  line=$(grep "np=" $log_file)

  if [[ $line =~ np=([0-9.]+)\ m=([0-9.]+)\ a=([0-9.]+)\ b=([0-9.]+) ]]; then
    rys_order=${BASH_REMATCH[1]} 
    chebyshev_order=${BASH_REMATCH[2]} 
    x_min=${BASH_REMATCH[3]}
    x_max=${BASH_REMATCH[4]}
  fi

  n=$((rys_order*(chebyshev_order+1)))

  output="rys_chebyshev_coeffs_${rys_order}_${i}.h"

  cat<<EOF | tee $output

  #ifndef RYS_CHEBYSHEV_COEFFS_${rys_order}_${i}_H
  #define RYS_CHEBYSHEV_COEFFS_${rys_order}_${i}_H
  const struct {
    const int rys_order = $rys_order;
    const int chebyshev_order = $chebyshev_order;
    const double x_min = ${x_min};
    const double x_max = ${x_max};
    const double roots_coefficients[$n] = {
$(egrep "^Coeffs: " $log_file | awk '{printf "%+38s, // %3d%3d%3d %4d\n", $5, $2, $3, $4, NR}')
  };
    const double weights_coefficients[$n] = {
$(egrep "^Coeffs: " $log_file | awk '{printf "%+38s, // %3d%3d%3d %4d\n", $6, $2, $3, $4, NR}')
  };
} rys_chebyshev_coeffs_${rys_order}_${i};
#endif /* RYS_CHEBYSHEV_COEFFS_${rys_order}_${i}_H */

EOF

done

for((i=0; i<10; i++)); do 
  echo \#include \"rys_chebyshev_coeffs_${rys_order}_$i.h\"
done

#In[6]:= In[6]:= np=5 m=25 a=0 b=10