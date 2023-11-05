////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Clean a row vector
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Clean a matrix.
/// This function "cleans" a matrix with a space, i.e., it adds suitable linear combinations
/// of the rows in @em sub to the rows of @em mat such that all pivot columns in @em mat are
/// zero. Both matrices must be over the same field and have the same number of columns.
/// The second matrix, @em sub, must be in echelon form. The cleaned matrix is reduced to
/// echelon form.
/// @return Rank of the cleaned matrix, or -1 on error.

int matClean(Matrix_t *mat, const Matrix_t *sub)
{
   int i;

   /* Check the arguments
      ------------------- */
   matValidate(MTX_HERE, mat);
   matValidate(MTX_HERE, sub);
   if ((mat->field != sub->field) || (mat->noc != sub->noc)) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
      return -1;
   }
   if (sub->pivotTable == NULL) {
      mtxAbort(MTX_HERE,"Subspace: %s",MTX_ERR_NOTECH);
      return -1;
   }

   /* Clean
      ----- */
   for (i = 0; i < mat->nor; ++i) {
      PTR m = matGetPtr(mat,i);
      ffCleanRow(m,sub->data,sub->nor,sub->noc, sub->pivotTable);
   }

   /* Reduce to echelon form
      ---------------------- */
   if (matEchelonize(mat) < 0) {
      return -1;
   }
   return mat->nor;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
