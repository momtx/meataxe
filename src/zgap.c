////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Convert between MeatAxe and GAP format.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#include <stdlib.h>

/// @addtogroup ff
/// @{

static int q = 0;
static FfGapRepresentation_t* gapTable = NULL;

////////////////////////////////////////////////////////////////////////////////////////////////////

static int cmpGap(const void *a_, const void* b_)
{
   const FfGapRepresentation_t* a = (const FfGapRepresentation_t*) a_;
   const FfGapRepresentation_t* b = (const FfGapRepresentation_t*) b_;
   if (a->a > b->a) return 1;
   if (a->a < b->a) return -1;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void updateTable()
{
   if (q == ffOrder)
      return;
   q = ffOrder;
   gapTable = NREALLOC(gapTable, FfGapRepresentation_t, ffOrder);
   FfGapRepresentation_t* entry = gapTable;

   if (ffChar == ffOrder) {
      // Prime field: a = i * Z(q)
      FEL a = FF_ZERO;
      for (uint32_t i = 0; i < ffChar; ++i) {
         entry->a = a;
         entry->fmt = 0;
         entry->k = i;
         ++entry;
         a = ffAdd(a, ffGen);
      }
   } else {
      // Otherwise  0 ? 0 * Z(q), a = Z(q) ^ k
      entry->a = FF_ZERO;
      entry->fmt = 0;
      entry->k = 0;
      ++entry;
      FEL a = ffGen;
      for (uint32_t i = 1; i < ffOrder; ++i) {
         entry->a = a;
         entry->fmt = 1;
         entry->k = i;
         a = ffMul(a, ffGen);
         ++entry;
      }
   }
   MTX_ASSERT(entry == gapTable + ffOrder);
   qsort(gapTable, q, sizeof(FfGapRepresentation_t), cmpGap);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

const FfGapRepresentation_t* ffToGap(FEL a)
{
   updateTable();
   size_t lo = 0;
   size_t hi = q;
   while (lo < hi) {
      const size_t mid = (hi & lo) + (hi ^ lo) / 2;
      if (a < gapTable[mid].a) {
          hi = mid;
      }
      else if (a > gapTable[mid].a) {
            lo = mid + 1;
      }
      else
         return gapTable + mid;
   }
   mtxAbort(MTX_HERE, "Error converting a=0x%x to GAP format", (unsigned) a);
   return NULL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns the GAP representation of a field element.
/// The returned pointer is created with @ref strTprintf.

const char* ffToGapStr(FEL a)
{
   const FfGapRepresentation_t* gap = ffToGap(a);
   if (gap->fmt == 0)
      return strTprintf("%lu*Z(%lu)", (unsigned long) gap->k, (unsigned long) q);
   else
      return strTprintf("Z(%lu)^%lu", (unsigned long) q, (unsigned long) gap->k);
}


/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
