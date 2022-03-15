#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstring>

#include "rys_chebyshev_coeffs.h"
#include "rys_utils.h"
#include "primitive_gaussian.h"
#include "rys_2e_integral.h"
#include "crys2.h"

#include "fort.h"
extern "C" void FORT(fmt)(double *fm, const double &x, const int &n);

int main()
{
  std::cout.precision(16);

  const PrimitiveGaussian s1(0, 0, 0, 0.1234, 0.7, 2.1, 1.2);
  const PrimitiveGaussian s2(0, 0, 0, 0.345, -1.0, 1.9, 1.0);
  const PrimitiveGaussian s3(0, 0, 0, 0.1245, 0.8, 2.02, 0.8);
  const PrimitiveGaussian p1(1, 0, 0, 0.2256, 0.8, 2.04, 1.0);
  const PrimitiveGaussian p2(0, 1, 0, 0.2256, 0.8, 2.04, 1.0);
  const PrimitiveGaussian p3(0, 0, 1, 0.2256, 0.8, 2.04, 1.0);

  std::cout << s1 << std::endl;
  std::cout << s2 << std::endl;
  std::cout << s3 << std::endl;
  std::cout << p1 << std::endl;
  std::cout << p2 << std::endl;
  std::cout << p3 << std::endl;

  Rys2EIntegral rys_1(s1, s2, p1, p2);
  std::cout << rys_1 << std::endl;

  std::cout << " " << coulomb_repulsion(s1, s2, p1, s3) << std::endl;
  std::cout << " " << coulomb_repulsion(s1, s2, p2, s3) << std::endl;
  std::cout << " " << coulomb_repulsion(s1, s2, p3, s3) << std::endl;

  const int m = 20;
  double *fm = new double [m];
  assert(fm);
  const double x = 1.0;

  FORT(fmt)(fm, x, m);

  std::cout << std::endl;
  for(int i = 0; i < m; i++)
    std::cout << " " << i << " " << fm[i] << std::endl;

  if(fm) { delete [] fm; fm = 0; }

  return 0;
}
