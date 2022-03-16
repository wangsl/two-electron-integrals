/* $Id: fns.h,v 1.1 2005/07/19 22:41:20 wangsl Exp $ */

#ifndef FNS_H
#define FNS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define UNDEF_VAL -1.0e16
#define SMALL_VAL 1.0e-10
#define UNDEF_INTVAL 314159265

#ifndef M_PI 
#define M_PI 3.1415926535897932385
#endif

inline double sq(double x) { return x*x; }
inline double cube(double x) { return x*x*x; }
inline double fourth(double x) { return sq(sq(x)); }
inline double fifth(double x) { return x*fourth(x); }
inline double sixth(double x) { return sq(cube(x)); }
inline double seventh(double x) { return x*sixth(x); }
inline double eighth(double x) { return sq(sq(sq(x))); }
inline double ninth(double x) { return cube(cube(x)); }
inline double tenth(double x) { return x*ninth(x); }
inline double eleventh(double x) { return x*tenth(x); }
inline double twelveth(double x) { return fourth(cube(x)); }
inline double thirteenth(double x) { return x*twelveth(x); }
inline double fourteenth(double x) { return sq(x)*twelveth(x); }
inline double fifteenth(double x) { return cube(fifth(x)); }
inline double sixteenth(double x) { return fourth(fourth(x)); }
inline double seventeenth(double x) { return x*sixteenth(x); }
inline double eighteenth(double x) { return sq(ninth(x)); }

inline void bound(double &x, double min, double max)
{
  assert(min < max);
  if (x < min)
    x = min;
  else if (x > max)
    x = max;
}

inline double factorial(int n)
{
  const double f[10] = {1,1,2,6,24,120,720,5040,40320,362880};
  assert(n > 0);
  if (n < 10)
    return f[n];
  double m;
  for (m = 1; n > 1; m *= n--)
    ;
  return m;
}

inline double double_factorial(int n)
{
  double d = 1.0;
  while (n > 1) {
    d *= n;
    n -= 2;
  }
  return d;
}

inline double combinations(int n, int k)
{
  assert(n >= k && k > 0);
  return factorial(n)/(factorial(k)*factorial(n-k));
}

inline int min(int x, int y) { return x < y ? x : y; }
inline int max(int x, int y) { return x > y ? x : y; }
inline double min(double x, double y) { return x < y ? x : y; }
inline double max(double x, double y) { return x > y ? x : y; }
inline double sim_round(double x) { return floor(x+0.5); }
inline double fractional_part(double x) { return x - sim_round(x); }
inline int is_almost_zero(double x) { return fabs(x) < SMALL_VAL; }
inline double zero_if_almost_zero(double x) { return is_almost_zero(x) ? 0 : x; }
inline int are_approximately_equal(double x, double y) { return is_almost_zero(x-y); }
inline int is_integer(double x) { return are_approximately_equal(x, sim_round(x)); }
inline int is_perfect_square(double x) { return is_integer(sqrt(x)); }
inline int is_perfect_cube(double x) { return is_integer(pow(x,1.0/3.0)); }
inline double periodic(double x, double l) { return x - l*sim_round(x/l); }
inline double ex_x(double x) { return exp(-x)/x; }
inline double ex2_x2(double x) { return ex_x(x*x); }

inline double legendre_polynomial(int n, double x)
{
  switch(n) {
  case 0:
    return 1;
  case 1:
    return x;
  case 2:
    return 0.5*(3*x*x - 1);
  case 3:
    return 0.5*(5*x*x*x - 3*x);
  default:
    assert(0);
  }
  return 0;
}

/* Given three points (x,y), returns the slope at the middle point */
inline double three_point_slope(double x1, double y1, double x2, double y2, double x3, double y3)
{
  return (sq(x3)*(y1-y2) - 2*x2*(x3*(y1-y2)+x1*(y2-y3)) + sq(x2)*(y1-y3) + sq(x1)*(y2-y3)) /
    ((x1-x2)*(x1-x3)*(x2-x3));
}

#endif
