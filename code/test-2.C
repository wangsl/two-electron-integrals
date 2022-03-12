#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstring>

#include "rys_chebyshev_coeffs.h"
#include "rys_utils.h"
#include "primitive_gaussian.h"

int main()
{
  std::cout.precision(14);

  const double x = 0.7;
  const double y = 2.1;
  const double z = 1.2;
  const double alpha = 0.1234;
  const int l = 0;
  const int m = 0;
  const int n = 0;

  const PrimitiveGaussian p(l, m, n, alpha, x, y, z);

  std::cout << p << std::endl;
  return 0;
}
