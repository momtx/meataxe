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
/// @a vecs. A new matrix is allocated and returned. @a noc is the number of
/// columns of the result, which must be a divisor of the number of columns
/// of @a vecs.
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
   if ((noc > vecs->Noc) || (vecs->Noc % noc != 0)) {
      mtxAbort(MTX_HERE,"noc=%d (vec:%d): %s",noc,vecs->Noc,MTX_ERR_BADARG);
      return NULL;
   }

   /* Convert the vector.
      ------------------- */
   result = matAlloc(vecs->Field,vecs->Noc / noc,noc);
   if (result == NULL) {
      return NULL;
   }
   for (i = 0; i < result->Nor; ++i) {
      if (matCopyRegion(result,i,0, vecs,n,i * noc,1,noc) != 0) {
         mtxAbort(MTX_HERE,"Copy failed");
	 return NULL;
      }
   }
   return result;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
