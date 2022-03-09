#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstring>
#include <math.h>

#include "rys_chebyshev_coeffs.h"
#include "gauss_hermite.h"

int main()
{
  RysChebyshev::setup_parameters();
  GaussHermite::setup_parameters();
  
  const int rys_order = 10;

  double *roots = new double [rys_order];
  assert(roots);
  double *weights = new double [rys_order];
  assert(weights);

  std::cout << std::endl;
  for(double x = 0.0; x < 120; x += 0.1) {
    RysChebyshev::calculate_rys_roots_and_weights(rys_order, x, roots, weights);

    for(int k = 0; k < rys_order; k++) 
      std::cout << " " 
                << k << " "   
                << std::fixed << std::setprecision(2) << x << " " 
                << std::scientific << std::setprecision(20) << roots[k] << " " << weights[k] << " "
                << std::endl;
    std::cout << std::endl;
  }

  if(roots) { delete [] roots; roots = 0; }
  if(weights) { delete [] weights; weights = 0; }

  return 0;
}
