
/* created at: 2022-03-12 18:44:09 */

#include <iostream>
using namespace std;
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "indent.h"
#include "primitive_gaussian.h"
#include "die.h"

ostream & operator <<(ostream &s, const PrimitiveGaussian &c)
{
  s << " {\n";
  IndentPush();
  c.write_fields(s);
  IndentPop();
  return s << Indent() << " }";
}

void PrimitiveGaussian::write_fields(ostream &s) const
{
  s << Indent() << "l " << l << "\n";
  s << Indent() << "m " << m << "\n";
  s << Indent() << "n " << n << "\n";
  s << Indent() << "alpha " << alpha << "\n";
  s << Indent() << "x " << x << "\n";
  s << Indent() << "y " << y << "\n";
  s << Indent() << "z " << z << "\n";
  s << Indent() << "norm " << norm << "\n";
}

