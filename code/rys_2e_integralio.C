
/* created at: 2022-03-16 23:19:43 */

#include <iostream>
using namespace std;
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "indent.h"
#include "rys_2e_integral.h"
#include "die.h"

ostream & operator <<(ostream &s, const Rys2EIntegral &c)
{
  s << " {\n";
  IndentPush();
  c.write_fields(s);
  IndentPop();
  return s << Indent() << " }";
}

void Rys2EIntegral::write_fields(ostream &s) const
{
  s << Indent() << "p0 " << p0 << "\n";
  s << Indent() << "p1 " << p1 << "\n";
  s << Indent() << "p2 " << p2 << "\n";
  s << Indent() << "p3 " << p3 << "\n";
}

