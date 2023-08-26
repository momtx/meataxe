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
   dest = matAlloc(src->Field,src->Noc,src->Nor);
   if (dest == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate result");
      return NULL;
   }
   d = dest->Data;
   for (i = 0; i < src->Noc; ++i) {
      int k;
      s = src->Data;
      for (k = 0; k < src->Nor; ++k) {
#if defined(MTX_DEBUG) && defined(PARANOID)
         FEL f;
         ffSetNoc(src->Noc);
         f = ffExtract(s,i);
         ffSetNoc(src->Nor);
         ffInsert(d,k,f);
#else
         ffInsert(d,k,ffExtract(s,i));
#endif
         s = (PTR)((char*) s + src->RowSize);
      }
      /*d = ffGetPtr(d,1,dest->Noc);*/
      d = (PTR)((char*) d + dest->RowSize);

   }
   return dest;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
