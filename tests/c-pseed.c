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
   long n = 0;
   for (int i = 0; i < 13; ++i) {
      PTR v = matGetPtr(vecs,i);
      n = MakeSeedVector(basis,n,v);
      ASSERT(n >= 0);	// no error
   }
   ASSERT(MakeSeedVector(basis,n,dummy) == -1); // no more seed vectors
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

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
