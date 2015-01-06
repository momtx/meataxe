
#ifndef CHECK_FUNCTION_TABLE

#include "meataxe.h"

Poly_t *RndPol(int fl, int mindeg, int maxdeg);

void TestPolAlloc(unsigned flags);
void TestPolAdd(unsigned flags);
void TestPolCompare(unsigned flags);
void TestPolGcd(unsigned flags);
void TestPolMul(unsigned flags);

#else

  { 601, "Polynomial allocation", TestPolAlloc },
  { 602, "Polynomial comparison", TestPolCompare },
  { 603, "Polynomial addition", TestPolAdd },
  { 604, "Polynomial Multiplication", TestPolMul },
  { 605, "Polynomial g.c.d.", TestPolGcd },

#endif
