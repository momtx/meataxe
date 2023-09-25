////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Multiplication of permutations
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup perm
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Multiply permutations.
/// This function multiplies @em dest from the right by @em src. Both
/// permutations must have the same degree.
/// @param dest Pointer to the first permutation.
/// @param src Pointer to the second permutation.
/// @return @em dest, or 0 on error.

Perm_t *permMul(Perm_t *dest, const Perm_t *src)
{
   permValidate(MTX_HERE, src);
   permValidate(MTX_HERE, dest);
   if (dest->degree != src->degree) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
   }

   uint32_t* d = dest->data;
   const uint32_t* const s = src->data;
   for (uint32_t i = dest->degree; i > 0; --i) {
      *d = s[*d];
      ++d;
   }
   return dest;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
