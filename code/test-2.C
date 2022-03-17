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

void libint_test();

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

  std::cout << std::endl
    << " Rys: " << coulomb_repulsion(s1, s2, p1, s3) << std::endl
    << " Rys: " << coulomb_repulsion(s1, s2, p2, s3) << std::endl
    << " Rys: " << coulomb_repulsion(s1, s2, p3, s3) << std::endl
    << std::endl;

  libint_test();
  std::cout << std::endl;

  for(int i = 0; i < 1; i++) {
    std::cout << std::endl << " " << i << std::endl;
    Rys2EIntegral rys2e_1(s1, s2, p1, s3);
    std::cout << " myRys: " << rys2e_1.coulomb_repulsion() << std::endl;

    Rys2EIntegral rys2e_2(s1, s2, p2, s3);
    std::cout << " myRys: " << rys2e_2.coulomb_repulsion() << std::endl;

    Rys2EIntegral rys2e_3(s1, s2, p3, s3);
    std::cout << " myRys: " << rys2e_3.coulomb_repulsion() << std::endl;
  }
  return 0;
}
