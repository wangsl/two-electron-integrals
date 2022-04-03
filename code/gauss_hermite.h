
#ifndef GAUSS_HERMITE_H
#define GAUSS_HERMITE_H

#include "gauss_hermite_rw.h"

#define MAX_GAUSS_HERMITE_PARAMETERS_LENGTH 20

struct GaussHermiteRootAndWeights {
  int order;
  double x_min;
  const double *roots;
  const double *weights;
};

namespace GaussHermite {
  static int has_setup_parameters = 0;
  static int parameters_length = MAX_GAUSS_HERMITE_PARAMETERS_LENGTH;
  static GaussHermiteRootAndWeights parameters[MAX_GAUSS_HERMITE_PARAMETERS_LENGTH];
  void setup_parameters();
  const GaussHermiteRootAndWeights *parameter(const int n_order);
};

#undef MAX_GAUSS_HERMITE_PARAMETERS_LENGTH

#endif /* GAUSS_HERMITE_H */