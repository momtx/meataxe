////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check functions for matrices.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdlib.h>
#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

static int testSeedVectors(Matrix_t* basis, Matrix_t* vecs, Matrix_t* cmp, PTR dummy)
{
   uint32_t n = 0;
   for (int i = 0; i < 13; ++i) {
      PTR v = matGetPtr(vecs,i);
      ASSERT_EQ_INT(svgMakeNext(v, &n, basis), 0);
   }
   ASSERT(svgMakeNext(dummy, &n, basis) == -1);       // no more seed vectors
   ASSERT(matCompare(vecs,cmp) == 0);           // seed vectors are as expected
   return 0;
}

TstResult SeedVectorGenerator()
{
   const int FIELD = 3;
   const int NOC = 3;
   SelectField(FIELD);
   Matrix_t* basis = matId(FIELD,NOC);
   Matrix_t* vecs = matAlloc(FIELD,13,NOC);
   Matrix_t* cmp = MkMat(13,NOC,
         1,0,0,  0,1,0,  1,1,0,  2,1,0,  0,0,1,
         1,0,1,  2,0,1,  0,1,1,  1,1,1,  2,1,1,
         0,2,1,  1,2,1,  2,2,1);
   PTR dummy = ffAlloc(1, NOC);

   int result = testSeedVectors(basis, vecs, cmp, dummy);
   matFree(vecs);
   matFree(cmp);
   matFree(basis);
   sysFree(dummy);
   return result;
}

TstResult SeedVectorGenerator_CheckLimits()
{
   {
      // ok, 2*17^7 - 1 < 2^32
      Matrix_t* basis = matId(17, 7);
      uint32_t vecno = 0;
      ASSERT_EQ_INT(svgMakeNext(NULL, &vecno, basis), 0);
      ASSERT_EQ_INT(vecno, 1U);
      matFree(basis);
   }

   {
      // failure, 2*17^8 - 1 >= 2^32
      Matrix_t* basis = matId(17, 8);
      uint32_t vecno = 0;
      ASSERT_ABORT(svgMakeNext(NULL, &vecno, basis));
      matFree(basis);
   }
   return 0;
}


// vim:fileencoding=utf8:sw=3:ts=8:et:cin
