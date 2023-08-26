////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Bit Strings, incidence relation
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

/// @addtogroup bs
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Bit string incidence relation.
/// This function returns 1 if and only if every bit which is set in @em a is also set in
/// @em b. Both bit strings must have the same size.
/// @return 1 if a⊆b, 0 if a⊈b, -1 on error.

int bsIsSub(const BitString_t *a, const BitString_t *b)
{
   register int i;
   register const unsigned long *ap, *bp;

   // check arguments
   bsValidate(MTX_HERE, a);
   bsValidate(MTX_HERE, b);
   if (a->Size != b->Size) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
      return -1;
   }

   // calculate result
   ap = (unsigned long const *) a->Data;
   bp = (unsigned long const *) b->Data;
   for (i = a->BufSize; i > 0; --i, ++ap, ++bp) {
      if ((*ap ^ (*ap & *bp)) != 0) {
         return 0;
      }
   }
   return 1;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
