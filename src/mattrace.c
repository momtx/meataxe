////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Trace of a matrix
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Trace of a Matrix.
/// This function calculates the sum of all diagonal elements of a matrix.
/// Note that the matrix need not be square.
/// @param mat Pointer to the matrix.
/// @return Trace of @a mat, (FEL)-1 on error.

FEL matTrace(const Matrix_t *mat)
{
   matValidate(MTX_HERE, mat);
   ffSetField(mat->field);
   FEL trace = FF_ZERO;
   PTR x = mat->data;
   const int maxi = mat->nor > mat->noc ? mat->noc : mat->nor;
   for (int i = 0; i < maxi; ++i) {
      trace = ffAdd(trace,ffExtract(x,i));
      ffStepPtr(&x, mat->noc);
   }
   return trace;
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
