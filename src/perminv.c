////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Inverse of a permutation
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup perm
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Inverse of a permutation
/// This function calulates the inverse of a permutation.
/// @param src Pointer to the permutation.
/// @return The inverse of @em src, or 0 on error.

Perm_t *permInverse(const Perm_t *src)
{
   register long i;
   register long *d, *s;
   Perm_t *inv;

   // Check arguments
   permValidate(MTX_HERE, src);

   // Allocate a new permutation
   inv = permAlloc(src->Degree);
   if (inv == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate result buffer");
      return NULL;
   }
   d = inv->Data;
   s = src->Data;

   // Invert
   for (i = src->Degree - 1, s = src->Data + i; i >= 0; --i, --s) {
      d[*s] = i;
   }

   return inv;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
