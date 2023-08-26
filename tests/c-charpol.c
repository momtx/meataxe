////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests for characteristic polynomial
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdlib.h>
#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

static int CheckPoly(Poly_t *p, int degree, ...)
{
   int i;
   va_list al;
   va_start(al,degree);
   ASSERT_EQ_INT(p->Degree,degree);
   for (i = 0; i <= degree; ++i) {
      ASSERT_EQ_INT(p->Data[i],ffFromInt(va_arg(al,int)));
   }
   va_end(al);
   polFree(p);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult CharacteristicPolynomial()
{
   int result = 0;

   SelectField(2);
   Matrix_t *a = MkMat(6,6,
             1,0,0,0,0,0,    0,1,1,0,0,0,    0,0,0,0,1,0,
             0,0,1,1,0,0,    0,0,0,0,0,1,    0,0,0,0,1,1);

   result |= CheckPoly(charPolFactor(a) ,1,1,1);
   result |= CheckPoly(charPolFactor(NULL),4,0,1,0,0,1);
   result |= CheckPoly(charPolFactor(NULL),1,1,1);
   ASSERT(charPolFactor(NULL) == NULL);
   matFree(a);
   return result;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
