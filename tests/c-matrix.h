#ifndef CHECK_FUNCTION_TABLE

#include "meataxe.h"

Matrix_t *RndMat(int fl, int nor, int noc);

void TestMatAlloc(unsigned flags);
void TestMatAddMul(unsigned flags);
void TestMatCompare(unsigned flags);
void TestMatCopy(unsigned flags);
void TestMatCut(unsigned flags);
void TestMatClean(unsigned flags);
void TestMatDup(unsigned flags);
void TestMatEchelonize(unsigned flags);
void TestMatId(unsigned flags);
void TestMatInv(unsigned flags);
void TestMatOrder(unsigned flags);
void TestNullSpace(unsigned flags);

#else

  { 311, "Matrix linear combination", TestMatAddMul },
  { 301, "Matrix order", TestMatOrder },
  { 302, "Matrix duplication", TestMatDup },
  { 303, "Matrix inversion", TestMatInv },
  { 304, "Matrix allocation", TestMatAlloc },
  { 305, "Gauss elimination", TestMatEchelonize },
  { 306, "Gauss elimination (2)", TestMatClean },
  { 307, "Matrix null-space", TestNullSpace },
  { 308, "Matrix comparison", TestMatCompare },
  { 309, "Matrix cutting", TestMatCut },
  { 310, "Matrix copying", TestMatCopy },
  { 312, "Identity matrix", TestMatId },

#endif
