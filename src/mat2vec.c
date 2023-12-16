////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Convert matrix to vector
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup tp
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Convert matrix to vector.
/// This function converts a matrix into a row vector by concatenating
/// the rows of the matrix. If @a mat is a r by c matrix, the resulting
/// vector has rc entries. Instead of allocating a new buffer for the
/// result, %matToVec() expects a pointer to a matrix, @a vecs, and puts
/// the vector into the @a n-th row of this matrix. Of course, @a vecs must
/// be over the smae field as @a mat, have rc columns and at least n+1 rows.
/// @see VectorToMatrix
/// @param mat Matrix to convert.
/// @param vecs Destination of the vector.
/// @param n Row number where the vector is stored.
/// @return 0 on success, -1 on error.

int MatrixToVector(const Matrix_t *mat, Matrix_t *vecs, int n)
{
   int i;
   matValidate(MTX_HERE, mat);
   matValidate(MTX_HERE, vecs);
   if ((mat->nor * mat->noc != vecs->noc)
       || (mat->field != vecs->field)) {
      mtxAbort(MTX_HERE,"mat and vecs: %s",MTX_ERR_INCOMPAT);
      return -1;
   }
   if ((n < 0) || (n >= vecs->nor)) {
      mtxAbort(MTX_HERE,"n=%d (nor=%d): %s",n,vecs->nor,MTX_ERR_BADARG);
      return -1;
   }

   // Convert the matrix.
   for (i = 0; i < mat->nor; ++i) {
      matCopyRegion(vecs,n,i * mat->noc, mat,i,0,1,mat->noc);
   }

   return 0;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
