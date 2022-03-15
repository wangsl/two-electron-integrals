
#ifndef CRYS2_H
#define CRYS2_H

#include "rys_chebyshev_coeffs.h"
#include "primitive_gaussian.h"

double coulomb_repulsion(const PrimitiveGaussian &p1, const PrimitiveGaussian &p2,
                         const PrimitiveGaussian &p3, const PrimitiveGaussian &p4);

#endif /* CRYS2_H */

