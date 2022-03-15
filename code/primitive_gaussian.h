#ifndef PRIMITIVE_GAUSSIAN_H
#define PRIMITIVE_GAUSSIAN_H

#include <iostream>

class PrimitiveGaussian 
{
public:
  PrimitiveGaussian(
    const int l_, const int m_, const int n_, 
    const double alpha_,
    const double x_, const double y_, const double z_, 
    double norm_ = -1.0);
  PrimitiveGaussian(const PrimitiveGaussian &);
  ~PrimitiveGaussian();

  const int l, m, n; // out
  const double alpha; // out
  double x, y, z; // out
  double norm; // out

  void update_coordinates(const double x_, const double y_, const double z_)
  { x = x_; y = y_; z = z_; }

private:
  /* IO */
  friend std::ostream & operator <<(std::ostream &s, const PrimitiveGaussian &c);
  void write_fields(std::ostream &s) const;
};

#endif /* PRIMITIVE_GAUSSIAN_H */