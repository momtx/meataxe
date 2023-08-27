////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check tensor product.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdlib.h>
#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestMatTensor2(int dim)
{
   int i;
   for (i = 0; i < 10; ++i) {
      int i1, i2, k1, k2;
      int nor1 = mtxRandomInt(dim);
      int nor2 = mtxRandomInt(dim);
      int noc1 = mtxRandomInt(dim);
      int noc2 = mtxRandomInt(dim);
      Matrix_t *m1 = RndMat(ffOrder,nor1,noc1);
      Matrix_t *m2 = RndMat(ffOrder,nor2,noc2);
      Matrix_t *m3 = matTensor(m1,m2);
      for (i1 = 0; i1 < nor1; ++i1) {
         PTR r1 = matGetPtr(m1,i1);
         for (i2 = 0; i2 < nor2; ++i2) {
            PTR r2 = matGetPtr(m2,i2);
            PTR r3 = matGetPtr(m3,i1 * nor2 + i2);
            for (k1 = 0; k1 < noc1; ++k1) {
               FEL f1 = ffExtract(r1,k1);
               for (k2 = 0; k2 < noc2; ++k2) {
                  FEL f2 = ffExtract(r2,k2);
                  FEL f3 = ffExtract(r3,k1 * noc2 + k2);
                  ASSERT_EQ_INT (ffMul(f1,f2), f3);
               }
            }
         }
      }
      matFree(m1);
      matFree(m2);
      matFree(m3);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Matrix_Tensor(int q)
{
    int result = 0;
    for (int dim = 1; result == 0 && dim < 50; dim += dim / 3 + 1) {
	result |= TestMatTensor2(dim);
    }
    return result;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
