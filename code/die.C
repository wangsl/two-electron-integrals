
/* $Id$ */

#include <unistd.h>
#include <cstdlib>
#include <iostream>
#include "die.h"

void die_at(const char *s, const char *file, int line)
{
  std::cout << " " << s << " at file: '" << file << "' line: " << line << "\n" << std::endl;
  std::cout.flush();
}

void die(const char *s)
{
  std::cout << s << "\n" << std::endl;
}



