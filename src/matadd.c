////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Add two matrices
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Sum of two matrices.
/// This function adds @em src to @em dest, overwriteing the previos value in @em dest.
/// The matrices must be over the same field and have the  same dimensions.
/// @return @em dest on success, 0 on error.

Matrix_t *matAdd(Matrix_t *dest, const Matrix_t *src)
{
   PTR dp, sp;
   register int n;

   /* Check arguments
      --------------- */
   matValidate(MTX_HERE,src);
   matValidate(MTX_HERE,dest);
   if ((dest->Field != src->Field) || (dest->Nor != src->Nor) ||
       (dest->Noc != src->Noc)) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
      return NULL;
   }

   /* Add <src> to <dest>
      ------------------- */
   dp = dest->Data;
   sp = src->Data;
   ffSetField(src->Field);
   for (n = src->Nor; n > 0; --n) {
      ffAddRow(dp,sp, src->Noc);
      ffStepPtr(&dp, src->Noc);
      ffStepPtr(&sp, src->Noc);
   }

   /* Delete the pivot table
      ---------------------- */
   mat_DeletePivotTable(dest);

   return dest;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
