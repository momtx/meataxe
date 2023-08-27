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
   ffSetField(mat->Field);
   FEL trace = FF_ZERO;
   PTR x = mat->Data;
   const int maxi = mat->Nor > mat->Noc ? mat->Noc : mat->Nor;
   for (int i = 0; i < maxi; ++i) {
      trace = ffAdd(trace,ffExtract(x,i));
      ffStepPtr(&x, mat->Noc);
   }
   return trace;
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
