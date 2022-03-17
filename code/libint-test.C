
#include <cstring>
#include <libint/libint.h>
#include "coord.h"
#include "fns.h"
#include "fort.h"

extern "C" void FORT(fmt)(double *fm, const double &x, const int &n);

inline double normalization_factor(const int nx, const int ny, const int nz, const double zeta)
{
  return pow(2.0*zeta/M_PI, 0.75) * pow(4.0*zeta, 0.5*(nx+ny+nz)) /
    sqrt(double_factorial(2*nx-1) * double_factorial(2*ny-1) * double_factorial(2*nx-1));
}

inline double renormalization_factor(
  const int A_nx, const int A_ny, const int A_nz, const double zeta_A,
  const int B_nx, const int B_ny, const int B_nz, const double zeta_B,
  const int C_nx, const int C_ny, const int C_nz, const double zeta_C,
  const int D_nx, const int D_ny, const int D_nz, const double zeta_D)
{
  double norm = 
    normalization_factor(A_nx, A_ny, A_nz, zeta_A) * 
    normalization_factor(B_nx, B_ny, B_nz, zeta_B) * 
    normalization_factor(C_nx, C_ny, C_nz, zeta_C) * 
    normalization_factor(D_nx, D_ny, D_nz, zeta_D);
  norm /= 
      normalization_factor(A_nx+A_ny+A_nz, 0, 0, zeta_A) * 
      normalization_factor(B_nx+B_ny+B_nz, 0, 0, zeta_B) * 
      normalization_factor(C_nx+C_ny+C_nz, 0, 0, zeta_C) * 
      normalization_factor(D_nx+D_ny+D_nz, 0, 0, zeta_D);
  return norm;
}

void auxiliary_integrals(
  double *fm, const int m, 
  const double rho, const double PQ2,
  const double S12, const double S34,
  const double C1, const double C2, const double C3, const double C4)
{
  FORT(fmt)(fm, rho*PQ2, m);
  const double prefactor = 2.0*sqrt(rho/M_PI)*S12*S34*C1*C2*C3*C4;
  for(int i = 0; i < m; i++)
    fm[i] *= prefactor;
}

inline double ss_overlap_integral(
  const double zeta_A, const Cartesian &A,
  const double zeta_B, const Cartesian &B)
{
  const double AB2 = (A-B).sq();
  const double zeta = zeta_A + zeta_B;
  const double sAB = pow(M_PI/zeta, 1.5)*exp(-zeta_A*zeta_B/zeta * AB2);
  return sAB;
}

inline int index(
  const int a, const int b, const int c, const int d,
  const int na, const int nb, const int nc, const int nd)
{ return ((a*nb+b)*nc+c)*nd+d; }

void libint_test()
{
  const Cartesian A(0.7, 2.1, 1.2);
  const Cartesian B(-1.0, 1.9, 1.0);
  const Cartesian C(0.8, 2.04, 1.0);
  const Cartesian D(0.8, 2.02, 0.8);

  const double zeta_A = 0.1234;
  const double zeta_B = 0.345;
  const double zeta_C = 0.2256;
  const double zeta_D = 0.1245;

  const double zeta = zeta_A + zeta_B;
  const double eta = zeta_C + zeta_D;
  const double rho = zeta*eta/(zeta+eta);
  const Cartesian P = (zeta_A*A + zeta_B*B)/zeta;
  const Cartesian Q = (zeta_C*C + zeta_D*D)/eta;
  const Cartesian W = (zeta*P + eta*Q)/(zeta+eta);

  const double S12 = ss_overlap_integral(zeta_A, A, zeta_B, B);
  const double S34 = ss_overlap_integral(zeta_C, C, zeta_D, D);

  const double C1 = normalization_factor(0, 0, 0, zeta_A);
  const double C2 = normalization_factor(0, 0, 0, zeta_B);
  const double C3 = normalization_factor(1, 0, 0, zeta_C);
  const double C4 = normalization_factor(0, 0, 0, zeta_D);

  const double PQ2 = (P-Q).sq();

  const int m = 2;
  double *fm = new double [m];
  assert(fm);
  auxiliary_integrals(fm, m, rho, PQ2, S12, S34, C1, C2, C3, C4);

  init_libint_base();

  Libint_t abcd;
  init_libint(&abcd, 1, 3);

  prim_data *data_array = new prim_data [3];
  assert(data_array);

  for(int k = 0; k < 3; k++) {
    prim_data &Data = data_array[k];
    Data.oo2z = 0.5/zeta;
    Data.oo2n = 0.5/eta;
    Data.poz = rho/zeta;
    Data.pon = rho/eta;
    Data.oo2zn = 0.5/(zeta+eta);
    Data.oo2p = 0.5/rho;
    memcpy(Data.F, fm, sizeof(double)*m);
    memcpy(Data.U[0], P-A, sizeof(Cartesian));
    memcpy(Data.U[2], Q-C, sizeof(Cartesian));
    memcpy(Data.U[4], W-P, sizeof(Cartesian));
    memcpy(Data.U[5], W-Q, sizeof(Cartesian));
  }

  abcd.PrimQuartet = data_array;
  Cartesian AB = A - B;
  Cartesian CD = C - D;
  memcpy(abcd.AB, AB, sizeof(Cartesian));
  memcpy(abcd.CD, CD, sizeof(Cartesian));

  double *a = build_eri[0][0][1][0](&abcd, 3);

  a[index(0, 0, 0, 0, 1, 1, 3, 1)] *= 
    renormalization_factor(0, 0, 0, zeta_A, 0, 0, 0, zeta_B,
                           1, 0, 0, zeta_C, 0, 0, 0, zeta_D);
  a[index(0, 0, 1, 0, 1, 1, 3, 1)] *= 
    renormalization_factor(0, 0, 0, zeta_A, 0, 0, 0, zeta_B,
                           0, 1, 0, zeta_C, 0, 0, 0, zeta_D);
  a[index(0, 0, 2, 0, 1, 1, 3, 1)] *= 
    renormalization_factor(0, 0, 0, zeta_A, 0, 0, 0, zeta_B,
                           0, 0, 1, zeta_C, 0, 0, 0, zeta_D);

  const double factor = 1.0/3.0;
  for(int i = 0; i < 3; i++)
    std::cout << " libint: " << a[index(0, 0, i, 0, 1, 1, 3, 1)]*factor << std::endl;

  // if(data_array) { delete [] data_array; data_array = 0; }
  if(fm) { delete [] fm; fm = 0; }
  free_libint(&abcd);
}