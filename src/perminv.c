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
   permValidate(MTX_HERE, src);

   Perm_t* inv = permAlloc(src->degree);
   uint32_t* d = inv->data;
   const uint32_t* s = src->data + src->degree;
   for (uint32_t i = src->degree; i > 0; ) {
      --i;
      --s;
      d[*s] = i;
   }

   return inv;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
