#include <iostream>
#include "gauss_hermite.h"

#define GAUSS_HERMITE(X) \
{ X.order, X.x_min, \
  (const double *) X.roots, \
  (const double *) X.weights \
}

void GaussHermite::setup_parameters()
{
  int i = -1;
  i++; GaussHermite::parameters[i] = GAUSS_HERMITE(gauss_hermite_2);
  i++; GaussHermite::parameters[i] = GAUSS_HERMITE(gauss_hermite_4);
  i++; GaussHermite::parameters[i] = GAUSS_HERMITE(gauss_hermite_6);
  i++; GaussHermite::parameters[i] = GAUSS_HERMITE(gauss_hermite_8);
  i++; GaussHermite::parameters[i] = GAUSS_HERMITE(gauss_hermite_10);
  i++; GaussHermite::parameters[i] = GAUSS_HERMITE(gauss_hermite_12);
  i++; GaussHermite::parameters[i] = GAUSS_HERMITE(gauss_hermite_14);
  i++; GaussHermite::parameters[i] = GAUSS_HERMITE(gauss_hermite_16);
  i++; GaussHermite::parameters[i] = GAUSS_HERMITE(gauss_hermite_18);
  i++; GaussHermite::parameters[i] = GAUSS_HERMITE(gauss_hermite_20);
  i++; GaussHermite::parameters[i] = GAUSS_HERMITE(gauss_hermite_22);
  i++; GaussHermite::parameters[i] = GAUSS_HERMITE(gauss_hermite_24);
  i++; GaussHermite::parameters[i] = GAUSS_HERMITE(gauss_hermite_26);
  i++; GaussHermite::parameters[i] = GAUSS_HERMITE(gauss_hermite_28);
  i++; GaussHermite::parameters[i] = GAUSS_HERMITE(gauss_hermite_30);
  i++; GaussHermite::parameters[i] = GAUSS_HERMITE(gauss_hermite_32);
  i++; GaussHermite::parameters[i] = GAUSS_HERMITE(gauss_hermite_34);
  i++; GaussHermite::parameters[i] = GAUSS_HERMITE(gauss_hermite_36);
  i++; GaussHermite::parameters[i] = GAUSS_HERMITE(gauss_hermite_38);
  i++; GaussHermite::parameters[i] = GAUSS_HERMITE(gauss_hermite_40);
}

const GaussHermiteRootAndWeights *GaussHermite::parameter(const int n_order)
  { return &parameters[n_order-1]; }

#undef GAUSS_HERMITE
