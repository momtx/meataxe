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

static int zmkpivot(PTR matrix, int nor, int noc, int *piv, int *ispiv)
{
   PTR x;
   int i, k;

   memset(ispiv,0,sizeof(int) * noc);

   /* Build the pivot table in <piv>.
      Keep track of assigned pivot columns in <ispiv>.
      ------------------------------------------------ */
   for (i = 0, x = matrix; i < nor && i < noc; ++i, ffStepPtr(&x, noc)) {
      FEL f;
      int newpiv = ffFindPivot(x,&f);
      if (ispiv[newpiv]) {
         mtxAbort(MTX_HERE, "%s", MTX_ERR_NOTECH);
         return -1;
      }
      piv[i] = newpiv;
      ispiv[newpiv] = 1;
   }

   /* Insert the non-pivot columns
      ---------------------------- */
   for (k = 0; k < noc; ++k) {
      if (!ispiv[k]) {
         piv[i++] = k;
      }
   }

   MTX_ASSERT(i == noc, -1);
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Reduce to echelon form.
/// This function builds the pivot table of a matrix. Unlike matEchelonize()
/// this function assumes that @a mat is already in echelon form.
/// @param mat Pointer to the matrix.
/// @return 0 on success, -1 on error.

int matPivotize(Matrix_t *mat)
{
   int *newtab;
   static int *is_pivot = NULL;
   static int maxnoc = -1;

   // check argument
   matValidate(MTX_HERE, mat);

   // (re-)allocate pivot table
   newtab = NREALLOC(mat->PivotTable,int,mat->Noc);
   if (newtab == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate pivot table (size %d)",mat->Noc);
      return -1;
   }
   mat->PivotTable = newtab;
   if (mat->Noc > maxnoc) {
      int *new_is_piv = NREALLOC(is_pivot,int,mat->Noc);
      if (new_is_piv == NULL) {
         mtxAbort(MTX_HERE,"Cannot allocate temporary table");
         return -1;
      }
      is_pivot = new_is_piv;
      maxnoc = mat->Noc;
   }

   // build the pivot table
   ffSetField(mat->Field);
   ffSetNoc(mat->Noc);
   return zmkpivot(mat->Data,mat->Nor,mat->Noc,mat->PivotTable,is_pivot);
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
