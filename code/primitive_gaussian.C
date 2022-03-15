
#include <cmath>
#include "primitive_gaussian.h"
#include "rys_utils.h"

PrimitiveGaussian::PrimitiveGaussian(
  const int l_, const int m_, const int n_, 
  const double alpha_,
  const double x_, const double y_, const double z_,
  double norm_) :
    l(l_), m(m_), n(n_), 
    alpha(alpha_), x(x_), y(y_), z(z_),
    norm(norm_)
  { 
    norm = pow(2*alpha/M_PI, 0.75) * pow(4*alpha, 0.5*(l+m+n)) /
          sqrt(double(double_factorial(2*l-1)) * double(double_factorial(2*m-1)) * double(double_factorial(2*n-1)));
  }

PrimitiveGaussian::PrimitiveGaussian(const PrimitiveGaussian &p) :
  l(p.l), m(p.m), n(p.n),
  alpha(p.alpha),
  x(p.x), y(p.y), z(p.z),
  norm(p.norm)
  { }

PrimitiveGaussian::~PrimitiveGaussian() 
{ }

std::ostream & operator <<(std::ostream &s, const PrimitiveGaussian &c)
{
  s << " { "
    << "(" << c.l << " " << c.m << " " << c.n << ") "
    << "(" << c.x << " " << c.y << " " << c.z << ") "
    << c.alpha << " " << c.norm 
    << " }"; 
  return s;
}