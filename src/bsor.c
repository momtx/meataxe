////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Bit Strings, "or" operation
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup bs
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Union of two bit strings.
/// This function computes the union of two bit strings, i.e., the bitwise logical "or".
/// The result is stored in @em dest and overwrites the previous value. Both bit
/// strings must have the same size.
/// @return 0 on success, -1 on error.

int bsOr(BitString_t *dest, const BitString_t *src)
{
   register int i;
   register unsigned long *dp;
   register const unsigned long *sp;

   // check arguments
   bsValidate(MTX_HERE, src);
   bsValidate(MTX_HERE, dest);
   if (dest->Size != src->Size) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
      return -1;
   }

   // OR operation
   dp = (unsigned long *) dest->Data;
   sp = (unsigned long const *) src->Data;
   for (i = src->BufSize; i > 0; --i) {
      *dp++ |= *sp++;
   }

   return 0;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
