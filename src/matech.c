////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Reduce a matrix to semi echelon form
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


////////////////////////////////////////////////////////////////////////////////////////////////////

static uint32_t zmkechelon(PTR matrix, uint32_t nor, uint32_t noc, uint32_t *piv, uint8_t *ispiv)
{
   uint32_t i;

   // Initialize the table
   for (i = 0; i < noc; ++i) {
      piv[i] = i;
      ispiv[i] = 0;
   }

   // Echelonize the matrix and build the pivot table in <piv>.
   // Keep track of assigned pivot columns in <ispiv>.
   PTR x = matrix;
   uint32_t rank = 0;
   PTR newrow = matrix;
   newrow = matrix;
   for (i = 0, x = matrix; i < nor && rank < noc; ++i, ffStepPtr(&x, noc)) {
      if (rank < i) {
         ffCopyRow(newrow,x, noc);
      }
      ffCleanRow(newrow,matrix,rank,noc,piv);
      FEL f;
      uint32_t newpiv = ffFindPivot(newrow,&f, noc);
      if (newpiv != MTX_NVAL) {
         piv[rank] = newpiv;
         ispiv[newpiv] = 1;
         ++rank;
         ffStepPtr(&newrow, noc);
      }
   }

   // Insert the non-pivot columns
   uint32_t j = rank;
   for (i = 0; i < noc; ++i) {
      if (!ispiv[i]) {
         piv[j++] = i;
      }
   }
   MTX_ASSERT(j == noc);

   return rank;
}


/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reduce to echelon form
/// This function performs a Gaussian elimination on the matrix @a mat.
/// After return, @a mat is in semi echelon form and a pivot table has been
/// attatched to the matrix. If the rank of @a mat was smaller than the number
/// of rows, some rows are removed during the process. This function can also
/// be used to rebuild the pivot table after the matrix has been modified.
/// @param mat Pointer to the matrix.
/// @return Rank of @a mat, or -1 on error.

int matEchelonize(Matrix_t *mat)
{
   matValidate(MTX_HERE, mat);

   // Re-allocate the pivot table. This is not really necessary, since
   // «Noc» should never change without releasing the pivot table, but
   // this would be a really nasty bug....
   mat->pivotTable = NREALLOC(mat->pivotTable,uint32_t,mat->noc);

   // Build the pivot table
   uint8_t *is_pivot = NALLOC(uint8_t, mat->noc);
   ffSetField(mat->field);
   const uint32_t rank = zmkechelon(mat->data,mat->nor,mat->noc,mat->pivotTable, is_pivot);
   sysFree(is_pivot);

   // If the rank is less than the number of rows, remove null rows
   if (rank != mat->nor) {
      mat->nor = rank;
      mat->data = (PTR) sysRealloc(mat->data,ffSize(rank, mat->noc));
   }

   return rank;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Nullity of a matrix.
/// This function calculates the dimension of the null-space of a matrix.
/// Unlike matNullity__() this function does not modify the matrix.
/// @param mat Pointer to the matrix.
/// @return Nullity of the matrix, or -1 on error.

int matNullity(const Matrix_t *mat)
{
   if (mat == NULL) {
      return -1;
   }
   return matNullity__(matDup(mat));
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Nullity of a matrix.
/// This function calculates the dimension of the null-space of a matrix
/// and deletes the matrix.
/// @param mat Pointer to the matrix.
/// @return Nullity of @em mat, or -1 on error.

int matNullity__(Matrix_t *mat)
{
   if (mat == NULL) {
      return -1;
   }
   if (matEchelonize(mat) < 0) {
      return -1;
   }
   int nul = mat->noc - mat->nor;
   matFree(mat);
   return nul;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
