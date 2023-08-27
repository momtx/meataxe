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
/// zero. Both matrices must be over the same field and have the same number of colums.
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
   if ((mat->Field != sub->Field) || (mat->Noc != sub->Noc)) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
      return -1;
   }
   if (sub->PivotTable == NULL) {
      mtxAbort(MTX_HERE,"Subspace: %s",MTX_ERR_NOTECH);
      return -1;
   }

   /* Clean
      ----- */
   for (i = 0; i < mat->Nor; ++i) {
      PTR m = matGetPtr(mat,i);
      ffCleanRow(m,sub->Data,sub->Nor,sub->Noc, sub->PivotTable);
   }

   /* Reduce to echelon form
      ---------------------- */
   if (matEchelonize(mat) < 0) {
      return -1;
   }
   return mat->Nor;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
