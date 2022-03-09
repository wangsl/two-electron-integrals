
#ifndef GAUSS_HERMITE_H
#define GAUSS_HERMITE_H

#include "gauss_hermite_rw.h"

struct GaussHermiteRootAndWeights {
  int order;
  double x_min;
  const double *roots;
  const double *weights;
};

namespace GaussHermite {
  static GaussHermiteRootAndWeights parameters[20];
  void setup_parameters();
  const GaussHermiteRootAndWeights *parameter(const int n_order);
};

#endif /* GAUSS_HERMITE_H */