////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Rebuild the pivot table of a matrix
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

static int zmkpivot(PTR matrix, uint32_t nor, uint32_t noc, uint32_t *piv, uint8_t *ispiv)
{
   // Extract the pivot columns to «piv».
   memset(ispiv,0,sizeof(uint8_t) * noc);
   PTR x = matrix;
   uint32_t i = 0;
   for (i = 0; i < nor && i < noc; ++i, ffStepPtr(&x, noc)) {
      FEL f;
      uint32_t newpiv = ffFindPivot(x,&f, noc);
      if (newpiv == MTX_NVAL || ispiv[newpiv])
         mtxAbort(MTX_HERE, "%s", MTX_ERR_NOTECH);
      piv[i] = newpiv;
      ispiv[newpiv] = 1;
   }

   // Append the non-pivot columns
   for (uint32_t k = 0; k < noc; ++k) {
      if (!ispiv[k]) {
         piv[i++] = k;
      }
   }

   MTX_ASSERT(i == noc);
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Create pivot table.
/// This function creates or updates the pivot table of a matrix. Unlike @c matEchelonize() this
/// function assumes that @a mat is already in echelon form. If it is not, @c matPivotize() fails
/// and aborts the program.

void matPivotize(Matrix_t *mat)
{
   matValidate(MTX_HERE, mat);

   mat->pivotTable = NREALLOC(mat->pivotTable,uint32_t,mat->noc);
   uint8_t *isPivot = NALLOC(uint8_t, mat->noc);
   ffSetField(mat->field);
   zmkpivot(mat->data,mat->nor,mat->noc,mat->pivotTable,isPivot);
   sysFree(isPivot);
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
