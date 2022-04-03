
#ifndef RYS_CHEBYSHEV_COEFFS_H
#define RYS_CHEBYSHEV_COEFFS_H

struct RysChebyshevCoeffs {
  int rys_order;
  int chebyshev_order;
  double x_min;
  double x_max;
  const double *roots_coefficients;
  const double *weights_coefficients;
};

#define MAX_RYS_CHEBYSHEV_COEFF_LENGTH 156

namespace RysChebyshev {
  //static int parameters_u = 1;
  static int has_setup_parameters = 0;
  static int parameters_length = MAX_RYS_CHEBYSHEV_COEFF_LENGTH;
  static RysChebyshevCoeffs parameters[MAX_RYS_CHEBYSHEV_COEFF_LENGTH];

  void setup_parameters();
  const RysChebyshevCoeffs *parameter(const int n_rys, const double x);
  void calculate_rys_roots_and_weights(const int n_rys, const double x, 
                                       double *roots, double *weights, 
                                       const int need_t = 0);
  void test();
  void initialize();
  static void setup_parameters_u();
};

#endif /* RYS_CHEBYSHEV_COEFFS_H */
