
#ifndef RYS_2E_INTEGRAL_H
#define RYS_2E_INTEGRAL_H

#include <iostream>
#include "primitive_gaussian.h"

#define MAX_RYS_ROOT 10

class Rys2EIntegral
{
public:
  Rys2EIntegral(const PrimitiveGaussian &p0_, const PrimitiveGaussian &p1_,
                const PrimitiveGaussian &p2_, const PrimitiveGaussian &p3_);

  const PrimitiveGaussian &p0; // out
  const PrimitiveGaussian &p1; // out 
  const PrimitiveGaussian &p2; // out 
  const PrimitiveGaussian &p3; // out

private:

  double B00, B1, B1p, C, Cp;
  double G[MAX_RYS_ROOT][MAX_RYS_ROOT];
  double roots[MAX_RYS_ROOT];
  double weights[MAX_RYS_ROOT];

  void recur_factors_gamess(const double t, const double A, const double B,
                            const double Px, const double Qx, const double xi, const double xk);

  void recur(const double t, 
    const int i, const int j, const int k, const int l,
    const double xi, const double xj, const double xk, const double xl,
    const double alphai, const double alphaj, const double alphak, const double alphal);

  double shift(const int i, const int j, const int k, const int l, const double xij, const double xkl);

  double int1d(const double t, 
    const int ix, const int jx, const int kx, const int lx,
	  const double xi, const double xj, const double xk, const double xl,
	  const double alphai, const double alphaj, const double alphak, const double alphal);

  friend std::ostream & operator <<(std::ostream &s, const Rys2EIntegral &c);
  void write_fields(std::ostream &s) const;
};

#endif /* RYS_2E_INTEGRAL_H */