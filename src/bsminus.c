////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Bit Strings, difference
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup bs
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Difference of two bit strings.
/// This function computes the (set theoretical) difference of two bit strings, i.e., a bit
/// in the result is set if and only if it is set in @em dest but not in @em src.
/// The result is stored in @em dest and overwrites the previous value.
/// The arguments must be bit strings of the same size.
/// @return 0 on success, -1 on error.

int bsMinus(BitString_t *dest, const BitString_t *src)
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

   // calculate result
   dp = (unsigned long *) dest->Data;
   sp = (unsigned long const *) src->Data;
   for (i = src->BufSize; i > 0; --i) {
      *dp++ &= ~*sp++;
   }

   return 0;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
