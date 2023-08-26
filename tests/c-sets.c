////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests for integer sets
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "testing.h"
#include "meataxe.h"

#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

#define NMAT 5

TstResult setAllocation()
{
   Set_t  *m[NMAT];
   int i;

   for (i = 0; i < NMAT; ++i) {
      m[i] = setAlloc();
   }

   for (i = 0; i < NMAT; ++i) {
      ASSERT_EQ_INT(setIsValid(m[i]), 1);
      ASSERT_EQ_INT(m[i]->Size, 0);
   }
   for (i = 0; i < NMAT; ++i) {
      ASSERT(setFree(m[i]) == 0);
   }

   for (i = 0; i < NMAT; ++i) {
      ASSERT(!setIsValid(m[i]));
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Set_BasicOperations()
{
   Set_t *s;
   long d[100];
   int i;

   memset(d,0,sizeof(d));
   for (i = 1; i <= 100; ++i) {
      int p;
      for (p = mtxRandomInt(100); d[p] != 0; p = (p + 1) % 100) {
      }
      d[p] = i;
   }

   s = setAlloc();
   for (i = 0; i < 100; ++i) {
      int k;
      setInsert(s,d[i]);
      ASSERT(s->Size == i + 1);
      for (k = 0; k <= i; ++k) {
         ASSERT(setContains(s,d[k]));
      }
      for (k = i + 1; k < 100; ++k) {
         ASSERT(!setContains(s,d[k]));
      }
   }
   setFree(s);
   return 0;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
