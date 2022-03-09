#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstring>
#include <math.h>

#include "rys_chebyshev_coeffs.h"
#include "gauss_hermite.h"

inline void chebyshev_t(const int n, const double x, double *t)
{
  assert(-1.0 <= x && x <= 1.0);

  if(n >= 1) t[0] = 1.0;
  if(n >= 2) t[1] = x;
  
  if(n <= 2) return;

  const double two_x = x + x; 
  for(int i = 2; i <= n; i++)   
    t[i] = two_x*t[i-1] - t[i-2];
}

void calculate_rys_roots_and_weights(const double x, const int rys_order)
{
  const GaussHermiteRootAndWeights *gh_parameger = GaussHermite::parameter(rys_order);
  if(x > gh_parameger->x_min) {
    const double x2 = 1.0/sqrt(x);
    for(int k = 0; k < rys_order; k++) {
      const double r_gh = x2*gh_parameger->roots[k];
      const double w_gh = x2*gh_parameger->weights[k];
      std::cout << k << " " 
              << std::fixed << std::setprecision(2) << x << " " 
              << std::scientific << std::setprecision(20) << r_gh << " " << w_gh 
              << std::endl;
    }
    return;
  }

  const RysChebyshevCoeffs *coeffs = RysChebyshev::parameter(rys_order, x);
  assert(coeffs && coeffs->rys_order == rys_order);

  const int &chebyshev_order = coeffs->chebyshev_order;
  const double &x_min = coeffs->x_min;
  const double &x_max = coeffs->x_max;
  
  double *T = new double [chebyshev_order+1];
  assert(T);
  chebyshev_t(chebyshev_order, (2*x-x_min-x_max)/(x_max-x_min), T);

  for(int k = 0; k < rys_order; k++) {
    const double *r_coeffs = coeffs->roots_coefficients + k*(chebyshev_order+1);
    const double *w_coeffs = coeffs->weights_coefficients + k*(chebyshev_order+1);
    double r_s = 0.0;
    double w_s = 0.0;
    for(int i = 0; i < chebyshev_order+1; i++) {
      r_s += r_coeffs[i]*T[i];
      w_s += w_coeffs[i]*T[i];
    }

    std::cout << k << " " 
              << std::fixed << std::setprecision(2) << x << " " 
              << std::scientific << std::setprecision(20) << r_s << " " << w_s << " "
              << std::endl;
  }

  if(T) { delete [] T; T = 0; }
}

int main()
{
  RysChebyshev::setup_parameters();
  GaussHermite::setup_parameters();

  const int rys_order = 5;

  for(double x=0.0; x<100; x+=0.05) {
    calculate_rys_roots_and_weights(x, rys_order);
    std::cout << std::endl;
  }

  return 0;
}
