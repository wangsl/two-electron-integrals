
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstring>

inline void chebyshev_t(const int n, const double x, double *t)
{
  assert(-1.0 <= x && x <= 1.0);

  if(n >= 1) t[0] = 1.0;
  if(n >= 2) t[1] = x;
  
  if(n <= 2) return;

  const double two_x = x + x; 
  for(int i = 2; i<=n; i++)   
    t[i] = two_x*t[i-1] - t[i-2];
}

int main()
{
  //const int np = std::cout.precision();
  //std::cout.precision(20);

  const double r[] = { 
    #include "rys-coeffs-5-0.txt"
  };

  const double w[] = { 
    #include "rys-weights-5-0.txt"
  };

  const int n_r = sizeof(r)/sizeof(double);
  std::cout << "nr = " << n_r << std::endl;

  const int n_w = sizeof(w)/sizeof(double);
  std::cout << "nw = " << n_w << std::endl;

  const int rys_order = 5;
  const int chebyshev_order = 25;
  const double x_min = 0.0;
  const double x_max = 10.0;

  double **r_coeffs = new double *[rys_order];
  assert(r_coeffs);
  double **w_coeffs = new double *[rys_order];
  assert(w_coeffs);
  for(int i = 0; i < rys_order; i++) {
    r_coeffs[i] = (double *) r + i*(chebyshev_order+1);
    w_coeffs[i] = (double *) w + i*(chebyshev_order+1);
  }

  const double x = 5.0;
  double *T = new double [chebyshev_order+1];
  memset(T, 0, sizeof(double)*(chebyshev_order+1));

  const double x_ = (2*x-x_min-x_max)/(x_max-x_min);
  chebyshev_t(chebyshev_order, x_, T);

  for(int k = 0; k < rys_order; k++) {
    double r_s = 0.0;
    double w_s = 0.0;
    for(int i = 0; i < chebyshev_order+1; i++) {
      r_s += r_coeffs[k][i]*T[i];
      w_s += w_coeffs[k][i]*T[i];
    }
    std::cout << k << " " 
              << std::fixed << std::setprecision(2) << x << " " 
              << std::scientific << std::setprecision(20) << r_s << " " << w_s 
              << std::endl;
  }

  if(r_coeffs) { delete [] r_coeffs; r_coeffs = 0; }
  if(w_coeffs) { delete [] w_coeffs; w_coeffs = 0; }
  if(T) { delete [] T; T = 0; }
  
  return 0;
}