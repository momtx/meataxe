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
      ASSERT(fpFree(p[i]) == 0);
      ASSERT(!fpIsValid(p[i]));
   }
   return 0;
}
