////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check functions for factored polynomials.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult FPoly_Allocation()
{
   const int NPOLY = 5;
   FPoly_t* p[NPOLY];
   int i;

   for (i = 0; i < NPOLY; ++i) {
      p[i] = fpAlloc(3);
      ASSERT_EQ_INT(p[i]->field, 3);
   }
   for (i = 0; i < NPOLY; ++i) {
      fpIsValid(p[i]);
   }
   for (i = 0; i < NPOLY; ++i) {
      fpFree(p[i]);
      ASSERT(!fpIsValid(p[i]));
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult FPoly_FailsOnWrongField()
{
   FPoly_t* fp = fpAlloc(2);
   Poly_t* p = polAlloc(3,0);
   ASSERT_ABORT(fpMulP(fp, p, 1));
   fpFree(fp);
   polFree(p);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult FPoly_OrderFactorsByDegree()
{
   FPoly_t* fp = fpAlloc(2);
   for (int deg = 5; deg >= 0; --deg) {
       Poly_t* p = polAlloc(2, deg);
       fpMulP(fp, p, 1);
       polFree(p);
   }
   ASSERT_EQ_INT(fp->nFactors, 6U);
   for (int i = 0; i <= 5; ++i) {
       ASSERT_EQ_INT(fp->factor[i]->degree, i);
   }

   fpFree(fp);
   return 0;
}
