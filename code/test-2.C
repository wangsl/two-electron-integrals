#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstring>

#include "rys_chebyshev_coeffs.h"
#include "rys_utils.h"

int main()
{
  //RysChebyshev::test();

  for(int n = 0; n < 30; n++)
    std::cout << " " << n << " " << fact(n) << std::endl;

  std::cout << binomial(12, 4) << std::endl;

  return 0;
}
