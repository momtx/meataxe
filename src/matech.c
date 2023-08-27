////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Reduce a matrix to semi echelon form
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


////////////////////////////////////////////////////////////////////////////////////////////////////

static int zmkechelon(PTR matrix, int nor, int noc, int *piv, int *ispiv)
{
   PTR x, newrow;
   int i, j, rank;

   /* Initialize the table
      -------------------- */
   for (i = 0; i < noc; ++i) {
      piv[i] = i;
      ispiv[i] = 0;
   }

   /* Echelonize the matrix and build the pivot table in <piv>.
      Keep track of assigned pivot columns in <ispiv>.
      --------------------------------------------------------- */
   rank = 0;
   newrow = matrix;
   for (i = 0, x = matrix; i < nor && rank < noc; ++i, ffStepPtr(&x, noc)) {
      int newpiv;
      FEL f;

      if (rank < i) {
         ffCopyRow(newrow,x, noc);
      }
      ffCleanRow(newrow,matrix,rank,noc,piv);
      newpiv = ffFindPivot(newrow,&f, noc);
      MTX_ASSERT(newpiv < noc, 0);
      if (newpiv >= 0) {
         piv[rank] = newpiv;
         ispiv[newpiv] = 1;
         ++rank;
         ffStepPtr(&newrow, noc);
      }
   }

   /* Insert the non-pivot columns
      ---------------------------- */
   j = rank;
   for (i = 0; i < noc; ++i) {
      if (!ispiv[i]) {
         piv[j++] = i;
      }
   }
   MTX_ASSERT(j == noc, 0);

   return rank;
}


/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Reduce to echelon form
/// This function performs a Gaussian elimination on the matrix |mat|. On
/// return, |mat| is in semi echelon form and a pivot table has been
/// attatched to the matrix. If the rank of |mat| was smaller than the number
/// of rows, some rows are removed during the process. This function can also
/// be used to rebuild the pivot table after the matrix has been modified.
/// @param mat Pointer to the matrix.
/// @return Rank of @em mat, or -1 on error.

int matEchelonize(Matrix_t *mat)
{
   int rank;
   int *newtab;
   static int *is_pivot = NULL;
   static int maxnoc = -1;

   /* Check the argument
      ------------------ */
   matValidate(MTX_HERE, mat);

   /* Re-allocate the pivot table. This is not really necessary, since
      |Noc| should never change without releasing the pivot table, but
      this would be a really nasty bug....
      ----------------------------------------------------------------- */
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

   /* Build the pivot table
      --------------------- */
   ffSetField(mat->Field);
   rank = zmkechelon(mat->Data,mat->Nor,mat->Noc,mat->PivotTable,is_pivot);

   /* If the rank is less than the number of rows, remove null rows
      ------------------------------------------------------------- */
   if (rank != mat->Nor) {
      mat->Nor = rank;
      mat->Data = (PTR) sysRealloc(mat->Data,ffSize(rank, mat->Noc));
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
   int nul = mat->Noc - mat->Nor;
   matFree(mat);
   return nul;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
