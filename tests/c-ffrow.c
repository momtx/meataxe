///////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check low-level row operations.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdlib.h>
#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestScalarProduct1(PTR a, PTR b, int size)
{
   int count;
   for (count = 0; count < 10; ++count) {
      int i;
      FEL expected = FF_ZERO;
      ffMulRow(a,FF_ZERO);
      ffMulRow(b,FF_ZERO);
      for (i = 0; i < size; ++i) {
         FEL f1 = FTab[mtxRandomInt(ffOrder)];
         FEL f2 = FTab[mtxRandomInt(ffOrder)];
         ffInsert(a,i,f1);
         ffInsert(b,i,f2);
         expected = ffAdd(expected,ffMul(f1,f2));
      }
      ASSERT_EQ_INT(ffScalarProduct(a,b), expected);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
 
TstResult ScalarProduct()
{
   int result = 0;
   for (int size = 0; result == 0 && size < 1000; size += size / 10 + 1) {
      PTR a, b;
      ffSetNoc(size);
      a = ffAlloc(1, size);
      b = ffAlloc(1, size);
      result |= TestScalarProduct1(a,b,size);
      sysFree(a);
      sysFree(b);
   }
   return result;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
