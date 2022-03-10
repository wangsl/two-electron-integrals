#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include "gauss_hermite.h"
#include "rys_chebyshev_coeffs.h"
#include "fort.h"

#define RYS_CHEBYSHEV(X) \
  { X.rys_order, X.chebyshev_order, \
    X.x_min, X.x_max, \
    (const double *) X.roots_coefficients, \
    (const double *) X.weights_coefficients \
  }

void RysChebyshev::setup_parameters()
{
  if(has_setup_parameters) return;

  std::cout << " RysChebyshev::setup_parameters" << std::endl;

  has_setup_parameters=1;

  int i=-1;
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_1_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_1_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_1_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_1_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_1_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_1_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_1_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_1_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_1_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_1_9);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_2_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_2_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_2_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_2_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_2_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_2_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_2_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_2_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_2_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_2_9);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_3_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_3_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_3_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_3_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_3_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_3_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_3_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_3_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_3_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_3_9);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_4_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_4_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_4_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_4_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_4_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_4_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_4_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_4_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_4_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_4_9);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_5_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_5_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_5_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_5_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_5_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_5_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_5_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_5_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_5_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_5_9);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_6_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_6_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_6_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_6_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_6_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_6_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_6_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_6_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_6_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_6_9);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_7_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_7_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_7_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_7_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_7_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_7_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_7_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_7_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_7_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_7_9);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_8_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_8_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_8_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_8_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_8_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_8_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_8_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_8_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_8_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_8_9);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_9_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_9_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_9_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_9_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_9_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_9_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_9_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_9_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_9_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_9_9);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_10_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_10_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_10_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_10_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_10_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_10_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_10_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_10_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_10_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_10_9);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_11_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_11_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_11_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_11_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_11_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_11_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_11_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_11_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_11_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_11_9);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_12_0);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_12_1);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_12_2);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_12_3);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_12_4);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_12_5);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_12_6);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_12_7);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_12_8);
  i++; RysChebyshev::parameters[i] = RYS_CHEBYSHEV(rys_chebyshev_coeffs_12_9);

  assert(i<= 120);
}

#undef RYS_CHEBYSHEV

const RysChebyshevCoeffs *RysChebyshev::parameter(const int n_rys, const double x)
{
  for(int i = 0; i < 120; i++) {
    const RysChebyshevCoeffs &p = parameters[i];
    if(n_rys == p.rys_order && p.x_min <= x && x <= p.x_max) 
      return &p;
  }
  return 0;
}

inline void chebyshev_t(const int n, const double x, double *t)
{
  assert(-1.0 <= x && x <= 1.0);

  if(n >= 1) t[0] = 1.0;
  if(n >= 2) t[1] = x;
  
  if(n <= 2) return;

  const double two_x = x + x; 
  for(int i = 2; i <= n; i++)   
    t[i] = two_x*t[i-1] - t[i-2];
}

void RysChebyshev::calculate_rys_roots_and_weights(
  const int rys_order, const double x, 
  double *roots, double *weights,
  const int need_u)
{
  RysChebyshev::setup_parameters();
  GaussHermite::setup_parameters();

  const GaussHermiteRootAndWeights *gh_parameger = GaussHermite::parameter(rys_order);
  if(x > gh_parameger->x_min) {
    const double x2 = 1.0/sqrt(x);
    for(int k = 0; k < rys_order; k++) {
      roots[k] = x2*gh_parameger->roots[k];
      weights[k] = x2*gh_parameger->weights[k];
    }
  } else {
    const RysChebyshevCoeffs *coeffs = RysChebyshev::parameter(rys_order, x);
    assert(coeffs && coeffs->rys_order == rys_order);

    const int &chebyshev_order = coeffs->chebyshev_order;
    const double &x_min = coeffs->x_min;
    const double &x_max = coeffs->x_max;
    
    double *T = new double [chebyshev_order+1];
    assert(T);
    chebyshev_t(chebyshev_order, (2*x-x_min-x_max)/(x_max-x_min), T);

    for(int k = 0; k < rys_order; k++) {
      const double *r_coeffs = coeffs->roots_coefficients + k*(chebyshev_order+1);
      const double *w_coeffs = coeffs->weights_coefficients + k*(chebyshev_order+1);
      double &r_sum = roots[k];
      double &w_sum = weights[k];
      r_sum = 0.0;
      w_sum = 0.0;
      for(int i = 0; i < chebyshev_order+1; i++) {
        r_sum += r_coeffs[i]*T[i];
        w_sum += w_coeffs[i]*T[i];
      }
    }
    if(T) { delete [] T; T = 0; }
  }

  if(need_u) {
    for(int k = 0; k < rys_order; k++) {
      double &t = roots[k];
      const double t_sq = t*t;
      roots[k] = t_sq/(1.0-t_sq);
    }
  }
}                                                    

// Reference: J. Comput. Phys. 21, 144 (1976) Table 1 for n=5
/*
    0.00   1  0.148874338981631216  0.295524224714752870E+00
    0.00   2  0.433395394129247213  0.269266719309996461E+00
    0.00   3  0.679409568299024436  0.219086362515982069E+00
    0.00   4  0.865063366688984536  0.149451349150580587E+00
    0.00   5  0.973906528517172076  0.666713443086881657E-01

    5.00   1  0.120616479067479812  0.224047067536327305E+00
    5.00   2  0.359993608978937729  0.123946126190236719E+00
    5.00   3  0.591831842524405793  0.389863670917377156E-01
    5.00   4  0.802341831901655866  0.763621205330385338E-02
    5.00   5  0.956727269752466358  0.109653673890797621E-02

   10.00   1  0.101239395075509955  0.182931707896803408E+00
   10.00   2  0.304887573873516171  0.809149721111935422E-01
   10.00   3  0.511824380911261478  0.152238928314807464E-01
   10.00   4  0.722257402679194249  0.114136525649064825E-02
   10.00   5  0.921039112895102874  0.354524106744085883E-04
*/

void RysChebyshev::test()
{
  const int rys_order = 5;

  double *roots = new double [rys_order];
  assert(roots);
  double *weights = new double [rys_order];
  assert(weights);

  std::cout << std::endl;
  for(double x = 0.0; x < 120; x += 0.1) {
    RysChebyshev::calculate_rys_roots_and_weights(rys_order, x, roots, weights, 0);

    for(int k = 0; k < rys_order; k++) {
      std::cout << " " 
                << k << " "   
                << std::fixed << std::setprecision(2) << x << " " 
                << std::scientific << std::setprecision(20) << roots[k] << " " << weights[k]
                << std::endl;
    }
    std::cout << std::endl;
  }

  if(roots) { delete [] roots; roots = 0; }
  if(weights) { delete [] weights; weights = 0; }
}

// Fortran version: CalculateRysRootsAndWeights
extern "C" {
  void FORT(calculaterysrootsandweights)(const int &rys_order, const double &x, 
                                        double *roots, double *weights, 
                                        const int &need_u)
  {
    RysChebyshev::calculate_rys_roots_and_weights(rys_order, x, roots, weights, need_u);
    std::cout.flush();
  }
}