////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check functions for integer matrices.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

static IntMatrix_t* randomMatrix(uint32_t nor, uint32_t noc)
{
   IntMatrix_t *a = imatAlloc(nor, noc);
   for (size_t i = 0; i < a->nor * a->noc; ++i) {
       a->data[i] = mtxRandomInt(65530) - 32765;
   }
   return a;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult IntMatrix_Allocation()
{
   #define NMAT 5
   static const uint32_t nor[NMAT] = { 0,0,1,1,9 };
   static const uint32_t noc[NMAT] = { 0,1,0,1,9 };
   IntMatrix_t *m[NMAT];
   int i;

   for (i = 0; i < NMAT; ++i) {
      m[i] = imatAlloc(nor[i],noc[i]);
   }
   for (i = 0; i < NMAT; ++i) {
      ASSERT(m[i]->nor == nor[i]);
      ASSERT(m[i]->noc == noc[i]);
      for (size_t k = 0; k < nor[i] * noc[i]; ++k) {
	  ASSERT(m[i]->data[k] == 0);
      }
   }
   for (i = 0; i < NMAT; ++i) {
      imatFree(m[i]);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult IntMatrix_ThrowsOnDoubleFree()
{
   IntMatrix_t *m = imatAlloc(20, 30);
   imatFree(m);
   ASSERT_ABORT(imatFree(m));
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult IntMatrix_Duplicate()
{
   IntMatrix_t *m = randomMatrix(20, 30);
   IntMatrix_t *copy = imatDup(m);
   ASSERT(copy != m);
   ASSERT(copy->nor = m->nor);
   ASSERT(copy->noc = m->noc);
   ASSERT(memcmp(copy->data, m->data, sizeof(uint32_t) * m->nor * m->noc) == 0);
   imatFree(copy);
   imatFree(m);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult IntMatrix_CompareFindsSingleDifferingMark()
{
   IntMatrix_t *a = randomMatrix(20, 30);
   IntMatrix_t *b = imatDup(a);
   MTX_ASSERT(imatCompare(a, b) == 0);

   for (size_t i = 0; i < a->nor * a->noc; ++i) {
       a->data[i] += 1;
       ASSERT(imatCompare(a, b) == 1);
       b->data[i] += 2;
       ASSERT(imatCompare(a, b) == -1);
       a->data[i] += 1;
       ASSERT(imatCompare(a, b) == 0);
   }
   imatFree(a);
   imatFree(b);
   return 0;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin

