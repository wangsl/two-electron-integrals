#ifndef PRIMITIVE_GAUSSIAN_H
#define PRIMITIVE_GAUSSIAN_H

#include <iostream>

class PrimitiveGaussian 
{
public:
  PrimitiveGaussian(
    const int &l_, const int &m_, const int &n_, 
    const double &alpha_,
    const double &x_, const double &y_, const double &z_);
  PrimitiveGaussian(const PrimitiveGaussian &);
  ~PrimitiveGaussian();

  const int &l, &m, &n; // out
  const double &alpha; // out
  const double &x, &y, &z; // out
  double norm; // out

private:
   /* IO */
  friend std::ostream & operator <<(std::ostream &s, const PrimitiveGaussian &c);
  void write_fields(std::ostream &s) const;

  // assignment operator
  PrimitiveGaussian & operator =(const PrimitiveGaussian &);
};

#endif /* PRIMITIVE_GAUSSIAN_H */