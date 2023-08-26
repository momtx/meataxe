////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Reduce a matrix to semi echelon form.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Duplicate a matrix
/// This function creates a copy of an existing matrix. The caller is
/// responsible for destroying the copy with matFree() when it is no
/// longer needed.
/// @return A copy of the source Matrix, or 0 on error.

Matrix_t *matDup(const Matrix_t *src)
{
   Matrix_t *m;

   matValidate(MTX_HERE, src);
   m = matAlloc(src->Field,src->Nor,src->Noc);
   if (m == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate matrix");
      return NULL;
   }
   memcpy(m->Data,src->Data,ffSize(src->Nor, src->Noc));
   return m;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
