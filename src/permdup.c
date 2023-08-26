////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Duplicate a permutation
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup perm
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Duplicate a permutation.
/// This function creates a copy of an existing permutation.
/// @param src Pointer to the permutation.
/// @return Pointer to a copy of @em src or 0 on error.

Perm_t *permDup(const Perm_t *src)
{
   permValidate(MTX_HERE, src);
   Perm_t *p = permAlloc(src->Degree);
   if (p == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate result");
      return NULL;
   }
   memcpy(p->Data,src->Data,(size_t) src->Degree * sizeof(long));
   return p;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
