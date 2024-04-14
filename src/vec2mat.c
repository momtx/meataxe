////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Convert vector to matrix.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
//   Local data


/// @addtogroup tp
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Convert vector to matrix.
/// This function converts a vector with m=rc entries into a r by c
/// matrix by filling the matrix from top to bottom and left to right with
/// the entries of the vector. The vector is taken as the n-th row of
/// @p vecs. A new matrix is allocated and returned. @p noc is the number of
/// columns of the result, which must be a divisor of the number of columns
/// of @p vecs.
/// @see MatrixToVector
/// @param vecs List of vectors.
/// @param n Number of the vector to convert.
/// @param noc Desired number of columns.
/// @return The result matrix or NULL on error.

Matrix_t *VectorToMatrix(Matrix_t *vecs, int n, int noc)
{
   int i;
   Matrix_t *result;

   /* Check arguments.
      ---------------- */
   matValidate(MTX_HERE, vecs);
   if ((noc > vecs->noc) || (vecs->noc % noc != 0)) {
      mtxAbort(MTX_HERE,"noc=%d (vec:%d): %s",noc,vecs->noc,MTX_ERR_BADARG);
      return NULL;
   }

   /* Convert the vector.
      ------------------- */
   result = matAlloc(vecs->field,vecs->noc / noc,noc);
   if (result == NULL) {
      return NULL;
   }
   for (i = 0; i < result->nor; ++i) {
      matCopyRegion(result,i,0, vecs,n,i * noc,1,noc);
   }
   return result;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
