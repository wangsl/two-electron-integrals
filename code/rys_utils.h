#ifndef RYS_UTILS_H
#define RYS_UTILS_H

#include <cassert>
#include <cmath>

inline double sq(const double x) { return x*x; }

inline long int factorial(const int n)
{
  const int n_max = 21;
  assert(n < n_max);
  const long int predfinds[] = 
    { 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 
      479001600, 6227020800, 87178291200, 1307674368000, 20922789888000, 
      355687428096000, 6402373705728000, 121645100408832000, 2432902008176640000 };
  return n <= 1 ? 1 : predfinds[n];
}

inline long int double_factorial(const int n)
{
  const int n_max = 34;
  assert(n < n_max);
  const long int predefinds[] = 
    { 1, 1, 2, 3, 8, 15, 48, 105, 384, 945, 3840, 10395, 46080, 135135, 
      645120, 2027025, 10321920, 34459425, 185794560, 654729075, 
      3715891200, 13749310575, 81749606400, 316234143225, 1961990553600, 
      7905853580625, 51011754393600, 213458046676875, 1428329123020800, 
      6190283353629375, 42849873690624000, 191898783962510625, 
      1371195958099968000, 6332659870762850625 };
  return n <= 1 ? 1 : predefinds[n];
}

inline double distance_sq(const double x1, const double y1, const double z1,
                          const double x2, const double y2, const double z2 )
{ return sq(x2-x1) + sq(y2-y1) + sq(z2-z1); }

inline double distance(const double x1, const double y1, const double z1,
                       const double x2, const double y2, const double z2 )
{ return sqrt(distance_sq(x1, y1, z1, x2, y2, z2)); }

inline long binomial(const int a, const int b) 
{ return factorial(a)/(factorial(b)*factorial(a-b)); }

inline double product_center_1d(const double alpha_a, const double x_a, 
                                const double alpha_b, const double x_b)
{ return (alpha_a*x_a + alpha_b*x_b)/(alpha_a + alpha_b); }

#endif /* RYS_UTILS_H */