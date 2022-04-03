
#include <cmath>
#include "rys_2e_integral.h"
#include "rys_utils.h"
#include "rys_chebyshev_coeffs.h"
#include "fns.h"

Rys2EIntegral::Rys2EIntegral(
  const PrimitiveGaussian &p0_, const PrimitiveGaussian &p1_,
  const PrimitiveGaussian &p2_, const PrimitiveGaussian &p3_) :
  p0(p0_), p1(p1_), p2(p2_), p3(p3_) 
{ 
  const int n = max(max(p0.l+p1.l, p0.m+p1.m), p0.n+p1.n) + 1;
  const int m = max(max(p2.l+p3.l, p2.m+p3.m), p2.n+p3.n) + 1;
  G = new double *[n]; assert(G);
  for(int i = 0; i < n; i++) {
    G[i] = new double [m];
    assert(G[i]);
  }
}

Rys2EIntegral::~Rys2EIntegral()
{
  const int n = max(max(p0.l+p1.l, p0.m+p1.m), p0.n+p1.n) + 1;
  const int m = max(max(p2.l+p3.l, p2.m+p3.m), p2.n+p3.n) + 1;
  for(int i = 0; i < n; i++)
    if(G[i]) { delete G[i]; G[i] = 0; }
  if(G) { delete [] G; G = 0; }
 }

void Rys2EIntegral::recur_factors_gamess(
  const double t, 
  const double A, const double B,
  const double Px, const double Qx, 
  const double xi, const double xk)
{
  /* Analogous versions taken from Gamess source code */
  const double fff = t/(A+B)/(1+t);
  B00 = 0.5*fff;
  B1 = 1/(2*A*(1+t)) + 0.5*fff;
  B1p = 1/(2*B*(1+t)) + 0.5*fff;
  C = (Px-xi)/(1+t) + (B*(Qx-xi)+A*(Px-xi))*fff;
  Cp = (Qx-xk)/(1+t) + (B*(Qx-xk)+A*(Px-xk))*fff;
}

void Rys2EIntegral::recur(
  const double t, 
  const int i, const int j, const int k, const int l,
  const double xi, const double xj, const double xk, const double xl,
	const double alphai, const double alphaj, const double alphak, const double alphal)
{
  /* Form G(n,m)=I(n,0,m,0) intermediate values for a Rys polynomial */
  const int n = i + j;
  const int m = k + l;
  const double A = alphai + alphaj;
  const double B = alphak + alphal;
  const double Px = (alphai*xi + alphaj*xj)/A;
  const double Qx = (alphak*xk + alphal*xl)/B;

  recur_factors_gamess(t, A, B, Px, Qx, xi, xk);

  /* ABD eq 11. */
  G[0][0] = M_PI*exp(-alphai*alphaj*pow(xi-xj,2)/(alphai+alphaj)
		     -alphak*alphal*pow(xk-xl,2)/(alphak+alphal))/sqrt(A*B);

  if(n > 0) G[1][0] = C*G[0][0];  /* ABD eq 15 */
  if(m > 0) G[0][1] = Cp*G[0][0]; /* ABD eq 16 */

  for(int a = 2; a < n+1; a++) G[a][0] = B1*(a-1)*G[a-2][0] + C*G[a-1][0];
  for(int b = 2; b < m+1; b++) G[0][b] = B1p*(b-1)*G[0][b-2] + Cp*G[0][b-1];

  if((m == 0) || (n == 0)) return;

  for(int a=1; a < n+1; a++){
    G[a][1] = a*B00*G[a-1][0] + Cp*G[a][0];
    for(int b = 2; b < m+1; b++)
      G[a][b] = B1p*(b-1)*G[a][b-2] + a*B00*G[a-1][b-1] + Cp*G[a][b-1];
  }
}

double Rys2EIntegral::shift(
  const int i, const int j, const int k, const int l, 
  const double xij, const double xkl)
{
/* Compute and  output I(i,j,k,l) from I(i+j,0,k+l,0) (G) */
  /*  xij = xi-xj, xkl = xk-xl */
  double ijkl = 0.0;
  for(int m = 0; m < l+1; m++){
    double ijm0 = 0.0;
    for(int n = 0; n < j+1; n++) /* I(i,j,m,0)<-I(n,0,m,0)  */
      ijm0 += binomial(j,n)*pow(xij, j-n)*G[n+i][m+k];
    ijkl += binomial(l,m)*pow(xkl, l-m)*ijm0; /* I(i,j,k,l)<-I(i,j,m,0) */
  } 
  return ijkl;
}

double Rys2EIntegral::int1d(
  const double t,
  const int ix, const int jx, const int kx, const int lx,
  const double xi, const double xj, const double xk, const double xl,
  const double alphai, const double alphaj, const double alphak, const double alphal)
{
  recur(t, ix, jx, kx, lx, xi, xj, xk, xl, alphai, alphaj, alphak, alphal);
  return shift(ix, jx, kx, lx, xi-xj, xk-xl);
}

inline double product_center_1D(double alpha_a, double x_a, double alpha_b, double x_b)
{ return (alpha_a*x_a + alpha_b*x_b)/(alpha_a+alpha_b); }

inline double dist2(double x1, double y1, double z1, double x2, double y2, double z2)
{ return sq(x1-x2) + sq(y1-y2) + sq(z1-z2); }

double Rys2EIntegral::coulomb_repulsion()
{
  const int norder = (
    p0.l + p0.m + p0.n + 
    p1.l + p1.m + p1.n + 
    p2.l + p2.m + p2.n + 
    p3.l + p3.m + p3.n)/2 + 1;

  const double A = p0.alpha + p1.alpha;
  const double B = p2.alpha + p3.alpha;
  const double rho = A*B/(A+B);

  const double xp = product_center_1D(p0.alpha, p0.x, p1.alpha, p1.x);
  const double yp = product_center_1D(p0.alpha, p0.y, p1.alpha, p1.y);
  const double zp = product_center_1D(p0.alpha, p0.z, p1.alpha, p1.z);
  const double xq = product_center_1D(p2.alpha, p2.x, p3.alpha, p3.x);
  const double yq = product_center_1D(p2.alpha, p2.y, p3.alpha, p3.y);
  const double zq = product_center_1D(p2.alpha, p2.z, p3.alpha, p3.z);
  const double rpq2 = dist2(xp, yp, zp, xq, yq, zq);

  const double X = rpq2*rho;

  double *roots = new double [norder]; assert(roots);
  double *weights = new double [norder]; assert(weights);
  RysChebyshev::calculate_rys_roots_and_weights(norder, X, roots, weights);

  double sum = 0.0;
  for(int i = 0; i < norder; i++) {
    const double &t = roots[i];
    const double Ix = int1d(t, p0.l, p1.l, p2.l, p3.l, p0.x, p1.x, p2.x, p3.x, p0.alpha, p1.alpha, p2.alpha, p3.alpha);
    const double Iy = int1d(t, p0.m, p1.m, p2.m, p3.m, p0.y, p1.y, p2.y, p3.y, p0.alpha, p1.alpha, p2.alpha, p3.alpha);
    const double Iz = int1d(t, p0.n, p1.n, p2.n, p3.n, p0.z, p1.z, p2.z, p3.z, p0.alpha, p1.alpha, p2.alpha, p3.alpha);
    sum += Ix*Iy*Iz*weights[i];
  }

  if(roots) { delete [] roots; roots = 0; }
  if(weights) { delete [] weights; weights = 0; }

  return 2*sqrt(rho/M_PI) * p0.norm * p1.norm * p2.norm * p3.norm * sum; /* ABD eq 5 & 9 */
}
