////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Transpose a matrix
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup mat
/// @{

/// Transpose a matrix.
/// @param src Pointer to the matrix.
/// @return Pointer to the transposed matrix or 0 on error.

Matrix_t *matTransposed(const Matrix_t *src)
{
   PTR s, d;
   long i;
   Matrix_t *dest;

   matValidate(MTX_HERE, src);
   dest = matAlloc(src->field,src->noc,src->nor);
   if (dest == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate result");
      return NULL;
   }
   d = dest->data;
   for (i = 0; i < src->noc; ++i) {
      int k;
      s = src->data;
      for (k = 0; k < src->nor; ++k) {
#if defined(MTX_DEBUG) && defined(PARANOID)
         FEL f;
         f = ffExtract(s,i);
         ffInsert(d,k,f);
#else
         ffInsert(d,k,ffExtract(s,i));
#endif
         ffStepPtr(&s, src->noc);
      }
      ffStepPtr(&d, dest->noc);

   }
   return dest;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
