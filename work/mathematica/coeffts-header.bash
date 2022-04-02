#!/bin/bash

function setup_one_file()
{
  local log_file="$1"
  local line=$(grep "np=" $log_file)

  if [[ $line =~ np=([0-9.]+)\ m=([0-9.]+)\ a=([0-9.]+)\ b=([0-9.]+) ]]; then
    local rys_order=${BASH_REMATCH[1]} 
    local chebyshev_order=${BASH_REMATCH[2]} 
    local x_min=${BASH_REMATCH[3]}
    local x_max=${BASH_REMATCH[4]}
  fi

  local n=$((rys_order*(chebyshev_order+1)))
  local i=$((x_min/10))

  cat<<EOF 

const struct {
  const int rys_order = $rys_order;
  const int chebyshev_order = $chebyshev_order;
  const double x_min = ${x_min};
  const double x_max = ${x_max};
  const double roots_coefficients[$n] = {
$(egrep "^Coeffs: " $log_file | awk '{printf "%+36s, // %3d%3d%3d %4d\n", $5, $2, $3, $4, NR}')
  };
  const double weights_coefficients[$n] = {
$(egrep "^Coeffs: " $log_file | awk '{printf "%+36s, // %3d%3d%3d %4d\n", $6, $2, $3, $4, NR}')
  };
} rys_chebyshev_coeffs_u_${rys_order}_${i};
EOF
}

for((i=0; i<=155; i++)); do
  slurm_log="slurm-17060946_${i}.out"
  if [ ! -e ${slurm_log} ]; then
    echo "${slurm_log} does not exist"
    exit 
  fi
  setup_one_file ${slurm_log}
done
