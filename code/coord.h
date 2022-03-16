/* $Id: coord.h,v 1.5 2005/07/22 00:32:12 wangsl Exp wangsl $ */

#ifndef COORD_H
#define COORD_H

#include <iomanip>
#include <cmath>
#include <iostream>
//using namespace std;

class Cartesian
{
public:
  double x, y, z;
  //Cartesian() { }
  Cartesian(double xx = 0.0, double yy = 0.0, double zz = 0.0) : x(xx), y(yy), z(zz) { }
  Cartesian(const Cartesian &r) : x(r.x), y(r.y), z(r.z) { }
  Cartesian & operator+=(const Cartesian &c)
  { x += c.x; y += c.y; z += c.z; return *this; }
  Cartesian & operator-=(const Cartesian &c)
  { x -= c.x; y -= c.y; z -= c.z; return *this; }
  Cartesian & operator*=(double r)
  { x *= r; y *= r; z *= r; return *this; }
  Cartesian & operator/=(double r)
  { x /= r; y /= r; z /= r; return *this; }
  Cartesian & zero()
  { x = 0; y = 0; z = 0; return *this; }
  Cartesian & set(double xx, double yy, double zz)
  { x = xx; y = yy; z = zz; return *this; }
  Cartesian operator-() const
  { return Cartesian(-x,-y,-z); }
  friend Cartesian operator+(const Cartesian &c1, const Cartesian &c2)
  { return Cartesian(c1.x + c2.x, c1.y + c2.y, c1.z + c2.z); }
  friend Cartesian operator-(const Cartesian &c1, const Cartesian &c2)
  { return Cartesian(c1.x - c2.x, c1.y - c2.y, c1.z - c2.z); }
  friend double operator*(const Cartesian &c1, const Cartesian &c2)
  { return c1.x * c2.x + c1.y * c2.y + c1.z * c2.z; }
  friend Cartesian operator*(double r, const Cartesian &c)
  { return Cartesian(r*c.x, r*c.y, r*c.z); }
  friend Cartesian operator*(const Cartesian &c, double r)
  { return r*c; }
  friend Cartesian operator/(const Cartesian &c, double r)
  { return Cartesian(c.x/r, c.y/r, c.z/r); }
  double sq() const
  { return x*x + y*y + z*z; }
  double magnitude() const
  { return sqrt(sq()); }
  Cartesian & scale_to_unit_magnitude()
  { return *this /= magnitude(); }
  Cartesian as_unit_vector() const
  { return *this / magnitude(); }
  double dot(const Cartesian &c) const
  { return x*c.x + y*c.y + z*c.z; }
  Cartesian cross(const Cartesian &c) const
  { return Cartesian(y*c.z - z*c.y, z*c.x - x*c.z, x*c.y - y*c.x); }
  Cartesian multiply(const Cartesian &c) const
  { return Cartesian(x*c.x, y*c.y, z*c.z); }
  Cartesian divide(const Cartesian &c) const
  { return Cartesian(x/c.x, y/c.y, z/c.z); }
  operator double *()
  { return (double *) this; }
  operator const double *() const
  { return (const double *) this; }
  Cartesian & apply(double (*f)(double))
  { x = f(x); y = f(y); z = f(z); return *this; }
  Cartesian map(double (*f)(double)) const
  { return Cartesian(f(x),f(y),f(z)); }

  /* By Shenglong Wang 04-23-2004 */
  double sum() const
    { return x+y+z; }
  
  int is_same_with(const Cartesian &r) 
  {
    return (*this-r).sq() < 1.0e-10 ? 1 : 0; 
  }
  
  static double distance(const Cartesian &c1, const Cartesian &c2)
  { return (c1 - c2).magnitude(); }
  static double angle(const Cartesian &a, const Cartesian &b, const Cartesian &c)
  { 
    double t = (a-b).as_unit_vector() * (c-b).as_unit_vector();
    if (t > 1)
      t = 1;
    else if (t < -1)
      t = -1;
    return acos(t);
  }
  static double dihedral(const Cartesian &a, const Cartesian &b, 
			 const Cartesian &c, const Cartesian &d)
  {
    Cartesian ab = a - b;
    Cartesian bc = b - c;
    Cartesian cd = c - d;
    Cartesian abcd = ab.cross(cd);
    Cartesian abbc = ab.cross(bc).as_unit_vector();
    Cartesian bccd = bc.cross(cd).as_unit_vector();
    double t = abbc*bccd;
    if (t > 1)
      t = 1;
    else if (t < -1)
      t = -1;
    return (abcd*bc > 0 ? 1.0 : -1.0) * acos(t);
  }
  
  // By Shenglong Wang
  Cartesian & operator *=(const Cartesian &c)
  { x *= c.x; y *= c.y; z *= c.z;  return *this; }
  Cartesian & operator /=(const Cartesian &c)
  { x /= c.x; y /= c.y; z /= c.z;  return *this; }
  Cartesian cartesian_multiply(const Cartesian &c)
  { return Cartesian(x*c.x, y*c.y, z*c.z); }
  Cartesian cartesian_divide(const Cartesian &c)
  { return Cartesian(x/c.x, y/c.y, z/c.z); }
};

/* Tensor elements are stored by column to be compatible with Fortran */
class Tensor
{
public:
  double xx, yx, zx, xy, yy, zy, xz, yz, zz;
  Tensor(double r = 0.0) // scalar -- create diagonal tensor
    : xx(r), yx(0), zx(0), xy(0), yy(r), zy(0), xz(0), yz(0), zz(r) { } 
  Tensor(const double &axx, const double &axy, const double &axz,
	 const double &ayx, const double &ayy, const double &ayz,
	 const double &azx, const double &azy, const double &azz) : 
    xx(axx), yx(ayx), zx(azx),
    xy(axy), yy(ayy), zy(azy),
    xz(axz), yz(ayz), zz(azz) { }
  Tensor(const Cartesian &a, const Cartesian &b) : // outer product
    xx(a.x*b.x), yx(a.y*b.x), zx(a.z*b.x),
    xy(a.x*b.y), yy(a.y*b.y), zy(a.z*b.y),
    xz(a.x*b.z), yz(a.y*b.z), zz(a.z*b.z)
  { }
  Cartesian &col1() { return (Cartesian &) xx; }
  Cartesian &col2() { return (Cartesian &) xy; }
  Cartesian &col3() { return (Cartesian &) xz; }
  const Cartesian &col1() const { return (const Cartesian &) xx; }
  const Cartesian &col2() const { return (const Cartesian &) xy; }
  const Cartesian &col3() const { return (const Cartesian &) xz; }
  Cartesian row1() const { return Cartesian(xx,xy,xz); }
  Cartesian row2() const { return Cartesian(yx,yy,yz); }
  Cartesian row3() const { return Cartesian(zx,zy,zz); }
  void zero() { xx = xy = xz = yx = yy = yz = zx = zy = zz = 0; }
  void set_row1(const Cartesian &c) { xx = c.x; xy = c.y; xz = c.z; }
  void set_row2(const Cartesian &c) { yx = c.x; yy = c.y; yz = c.z; }
  void set_row3(const Cartesian &c) { zx = c.x; zy = c.y; zz = c.z; }
  operator double *()
  { return (double *) this; }
  operator const double *() const
  { return (const double *) this; }
  Tensor & operator+=(const Tensor d)
  {
    xx += d.xx; xy += d.xy; xz += d.xz;
    yx += d.yx; yy += d.yy; yz += d.yz;
    zx += d.zx; zy += d.zy; zz += d.zz;
    return *this;
  }
  Tensor & operator-=(const Tensor d)
  {
    xx -= d.xx; xy -= d.xy; xz -= d.xz;
    yx -= d.yx; yy -= d.yy; yz -= d.yz;
    zx -= d.zx; zy -= d.zy; zz -= d.zz;
    return *this;
  }
  Tensor & operator*=(double r)
  {
    xx *= r; xy *= r; xz *= r;
    yx *= r; yy *= r; yz *= r;
    zx *= r; zy *= r; zz *= r;
    return *this;
  }
  double trace() const { return xx + yy + zz; }
  Tensor transpose() const { return Tensor(xx,yx,zx,xy,yy,zy,xz,yz,zz); }
  double determinant() const { return col1()*col2().cross(col3()); }
  Cartesian diagonal() const { return Cartesian(xx,yy,zz); }
  Tensor inverse() const
  {
    Tensor t;
    t.set_row1(col2().cross(col3()));
    t.set_row2(col3().cross(col1()));
    t.set_row3(col1().cross(col2()));
    t *= 1/determinant();
    return t;
  }
  Tensor operator-() const { return Tensor(-xx,-xy,-xz,-yx,-yy,-yz,-zx,-zy,-zz); }
  Tensor & apply(double (*f)(double))
  {
    xx = f(xx); xy = f(xy); xz = f(xz);
    yx = f(yx); yy = f(yy); yz = f(yz);
    zx = f(zx); zy = f(zy); zz = f(zz);
    return *this;
  }
  Tensor map(double (*f)(double)) const
  { 
    return Tensor(f(xx),f(xy),f(xz),
		  f(yx),f(yy),f(yz),
		  f(zx),f(zy),f(zz)); 
  }
  friend Cartesian operator*(const Tensor &t, const Cartesian &c)
  {
    return Cartesian(t.xx*c.x + t.xy*c.y + t.xz*c.z,
		     t.yx*c.x + t.yy*c.y + t.yz*c.z,
		     t.zx*c.x + t.zy*c.y + t.zz*c.z);
  }
  friend Cartesian operator*(const Cartesian &c, const Tensor &t)
  {
    return Cartesian(c.x*t.xx + c.y*t.yx + c.z*t.zx,
		     c.x*t.xy + c.y*t.yy + c.z*t.zy,
		     c.x*t.xz + c.y*t.yz + c.z*t.zz);
  }
  friend Tensor operator*(double r, const Tensor &t)
  {
    return Tensor(r*t.xx, r*t.xy, r*t.xz,
		  r*t.yx, r*t.yy, r*t.yz,
		  r*t.zx, r*t.zy, r*t.zz);
  }
  friend Tensor operator/(const Tensor &t, double r)
  {
    return Tensor(t.xx/r, t.xy/r, t.xz/r,
		  t.yx/r, t.yy/r, t.yz/r,
		  t.zx/r, t.zy/r, t.zz/r);
  }
  friend Tensor operator*(const Tensor &t, double r)
  { return r*t; }
  friend Tensor operator+(const Tensor &t1, const Tensor &t2)
  {
    return Tensor(t1.xx+t2.xx, t1.xy+t2.xy, t1.xz+t2.xz,
		  t1.yx+t2.yx, t1.yy+t2.yy, t1.yz+t2.yz,
		  t1.zx+t2.zx, t1.zy+t2.zy, t1.zz+t2.zz);
  }
  friend Tensor operator-(const Tensor &t1, const Tensor &t2)
  {
    return Tensor(t1.xx-t2.xx, t1.xy-t2.xy, t1.xz-t2.xz,
		  t1.yx-t2.yx, t1.yy-t2.yy, t1.yz-t2.yz,
		  t1.zx-t2.zx, t1.zy-t2.zy, t1.zz-t2.zz);
  }
  friend Tensor operator*(const Tensor &d1, const Tensor &d2)
  {
    return Tensor(d1.xx*d2.xx + d1.xy*d2.yx + d1.xz*d2.zx,
		  d1.xx*d2.xy + d1.xy*d2.yy + d1.xz*d2.zy, 
		  d1.xx*d2.xz + d1.xy*d2.yz + d1.xz*d2.zz,
		  d1.yx*d2.xx + d1.yy*d2.yx + d1.yz*d2.zx, 
		  d1.yx*d2.xy + d1.yy*d2.yy + d1.yz*d2.zy, 
		  d1.yx*d2.xz + d1.yy*d2.yz + d1.yz*d2.zz,
		  d1.zx*d2.xx + d1.zy*d2.yx + d1.zz*d2.zx, 
		  d1.zx*d2.xy + d1.zy*d2.yy + d1.zz*d2.zy,
		  d1.zx*d2.xz + d1.zy*d2.yz + d1.zz*d2.zz);
  }
  static Cartesian z_location(const Cartesian &a, const Cartesian &b, const Cartesian &c,
			      double r, double theta, double phi)
  {
    Tensor j;
    j.col2() = (c-a).cross(b-a).as_unit_vector();
    j.col1() = (b-a).cross(j.col2()).as_unit_vector();
    j.col3() = j.col1().cross(j.col2()).as_unit_vector();
    return a + j * Cartesian(r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), -r*cos(theta));
  }

  double contract_with(const Tensor &a) {
    double s = 
      xx*a.xx + yx*a.yx + zx*a.zx +
      xy*a.xy + yy*a.yy + zy*a.zy +
      xz*a.xz + yz*a.yz + zz*a.zz;
    return s;
  }
  
  static Tensor dipole_dipole_tensor(const Cartesian &R)
  {
    const double r2 = R.sq();
    const double r3 = r2*sqrt(r2);
    const double r31 = 1.0/r3;
    const double r5 = r2*r3;
    Tensor T(R,R);
    T *= 3.0/r5;
    T.xx -= r31;
    T.yy -= r31;
    T.zz -= r31;
    return T;
  }

};

inline std::ostream & operator<<(std::ostream &s, const Cartesian &c) 
{ 
  /*
  char *buf = new char [256];
  sprintf(buf, "%20.10f%20.10f%20.10f", c.x, c.y, c.z);
  s << buf;
  delete [] buf;
  return s;
  */
  return s << c.x << " " << c.y << " " << c.z; 
}

inline std::istream & operator>>(std::istream &s, Cartesian &c)
{ 
  return s >> c.x >> c.y >> c.z; 
}

inline std::ostream & operator<<(std::ostream &s, const Tensor &c)
{ 
  /*
  char *buf = new char [256];
  s << "\n";
  sprintf(buf, "%20.10f%20.10f%20.10f\n", c.xx, c.xy, c.xz);
  s << buf;
  sprintf(buf, "%20.10f%20.10f%20.10f\n", c.yx, c.yy, c.yz);
  s << buf;
  sprintf(buf, "%20.10f%20.10f%20.10f", c.zx, c.zy, c.zz);
  s << buf;
  delete [] buf;
  return s;
  */
  s << c.xx << " " << c.xy << " " << c.xz << "\n"
    << c.yx << " " << c.yy << " " << c.yz << "\n"
    << c.zx << " " << c.zy << " " << c.zz << "\n";
  //<< "\n";
  return s;
}

inline std::istream & operator>>(std::istream &s, Tensor &c)
{ 
  return s >> c.xx >> c.xy >> c.xz
	   >> c.yx >> c.yy >> c.yz
	   >> c.zx >> c.zy >> c.zz;
}

#endif
