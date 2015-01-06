
#ifndef CHECK_FUNCTION_TABLE

#include "meataxe.h"

Perm_t *RndPerm(int degree);

void TestPermAlloc(unsigned flags);
void TestPermMul(unsigned flags);
void TestPermOrder(unsigned flags);
void TestPermPwr(unsigned flags);
void TestPermInv(unsigned flags);

#else

  { 501, "Permutation allocation", TestPermAlloc },
  { 502, "Permutation order", TestPermOrder },
  { 503, "Permutation multiplication", TestPermMul },
  { 504, "Permutation power", TestPermPwr },
  { 505, "Permutation inverse", TestPermInv },

#endif
