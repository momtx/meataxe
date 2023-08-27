///////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check low-level row operations.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdlib.h>
#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
 
TstResult ScalarProduct_RandomValues(int q)
{
   int result = 0;
   for (int noc = 0; result == 0 && noc < 30; ++noc) {
      PTR a = ffAlloc(1, noc);
      PTR b = ffAlloc(1, noc);
      ffMulRow(a,FF_ZERO, noc);
      ffMulRow(b,FF_ZERO, noc);
      FEL expected = FF_ZERO;
      for (int i = 0; i < noc; ++i) {
         FEL f1 = FTab[mtxRandomInt(ffOrder)];
         FEL f2 = FTab[mtxRandomInt(ffOrder)];
         ffInsert(a,i,f1);
         ffInsert(b,i,f2);
         expected = ffAdd(expected,ffMul(f1,f2));
      }
      FEL result = ffScalarProduct(a,b,noc);
      sysFree(a);
      sysFree(b);

      ASSERT_EQ_INT(result, expected);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
 
TstResult ScalarProduct_WorksForNocEqualsZero(int q)
{
   PTR a = ffAlloc(1, 0);
   PTR b = ffAlloc(1, 0);
   FEL result = ffScalarProduct(a, b, 0);
   sysFree(a);
   sysFree(b);
   ASSERT_EQ_INT(result, FF_ZERO);
   return 0;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
