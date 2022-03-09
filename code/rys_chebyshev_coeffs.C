#include <cassert>
#include "rys_chebyshev_coeffs.h"

#define RYS_CHEBYSHEV(X) \
  { X.rys_order, X.chebyshev_order, \
    X.x_min, X.x_max, \
    (const double *) X.roots_coefficients, \
    (const double *) X.weights_coefficients \
  }

void RysChebyshev::setup_parameters()
{
  int i=-1;
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_1_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_1_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_1_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_1_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_1_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_1_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_1_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_1_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_1_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_1_9);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_2_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_2_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_2_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_2_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_2_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_2_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_2_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_2_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_2_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_2_9);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_3_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_3_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_3_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_3_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_3_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_3_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_3_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_3_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_3_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_3_9);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_4_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_4_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_4_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_4_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_4_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_4_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_4_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_4_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_4_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_4_9);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_5_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_5_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_5_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_5_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_5_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_5_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_5_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_5_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_5_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_5_9);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_6_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_6_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_6_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_6_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_6_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_6_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_6_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_6_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_6_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_6_9);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_7_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_7_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_7_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_7_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_7_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_7_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_7_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_7_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_7_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_7_9);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_8_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_8_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_8_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_8_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_8_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_8_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_8_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_8_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_8_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_8_9);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_9_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_9_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_9_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_9_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_9_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_9_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_9_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_9_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_9_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_9_9);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_10_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_10_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_10_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_10_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_10_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_10_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_10_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_10_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_10_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_10_9);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_11_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_11_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_11_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_11_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_11_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_11_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_11_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_11_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_11_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_11_9);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_12_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_12_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_12_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_12_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_12_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_12_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_12_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_12_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_12_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_12_9);

  assert(i<= 120);
}

const RysChebyshevCoeffs *RysChebyshev::parameter(const int n_rys, const double x)
{
  for(int i = 0; i < 120; i++) {
    const RysChebyshevCoeffs &p = parameters[i];
    if(n_rys == p.rys_order && p.x_min <= x && x <= p.x_max) 
      return &p;
  }
  return 0;
}

#undef RYS_CHEBYSHEV