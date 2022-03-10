#ifndef RYS_UTILS_H
#define RYS_UTILS_H

#include <cmath>

inline double sq(const double x) { return x*x; }

inline long int fact(const int n)
{
  const int n_max = 21;
  const long int predfined_facts[] = 
    { 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 
      479001600, 6227020800, 87178291200, 1307674368000, 20922789888000, 
      355687428096000, 6402373705728000, 121645100408832000, 2432902008176640000 };
  if(n <= 1) return 1;
  return n < n_max ? predfined_facts[n] : -1;
}

inline double distance_sq(const double x1, const double y1, const double z1,
                          const double x2, const double y2, const double z2 )
{ return sq(x2-x1) + sq(y2-y1) + sq(z2-z1); }

inline double distance(const double x1, const double y1, const double z1,
                       const double x2, const double y2, const double z2 )
{ return sqrt(distance_sq(x1, y1, z1, x2, y2, z2)); }

inline long binomial(const int a, const int b) { return fact(a)/(fact(b)*fact(a-b)); }

inline double product_center_1d(const double alpha_a, const double x_a, 
                                const double alpha_b, const double x_b)
{ return (alpha_a*x_a + alpha_b*x_b)/(alpha_a + alpha_b); }

#endif /* RYS_UTILS_H */