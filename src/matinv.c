////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Matrix inversion
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


////////////////////////////////////////////////////////////////////////////////////////////////////

static void zmatinv(PTR mat, PTR result, int noc)
{
   PTR xj1, xj2, xk1, xk2;
   FEL f1 = FF_ZERO, f2;
   long j, k;

   // initialize result with identity matrix
   for (j = 0, xj1 = result; j < noc; ++j, ffStepPtr(&xj1, noc)) {
      ffMulRow(xj1,FF_ZERO, noc);
      ffInsert(xj1,j,FF_ONE);
   }

   // matrix inversion
   xj1 = mat;
   xj2 = result;
   for (j = 0; j < noc; ++j) {

      for (xk1 = xj1, k = j;
           k < noc && (f1 = ffExtract(xk1,j)) == FF_ZERO;
           ++k, ffStepPtr(&xk1, noc)) {
      }
      if (f1 == FF_ZERO) {
         mtxAbort(MTX_HERE,"%s",MTX_ERR_DIV0);
      }
      if (k > j) {
         ffSwapRows(xk1,xj1,noc);
         // xk2 = ffGetPtr(xj2,k-j,noc)
         xk2 = (PTR)((char *)xj2 + (k - j) * ffRowSize(noc));
         ffSwapRows(xk2,xj2,noc);
      }
      f2 = ffInv(f1);
      ffMulRow(xj1,f2,noc);
      ffMulRow(xj2,f2,noc);
      xk1 = mat;
      xk2 = result;
      for (k = 0; k < noc; ++k) {
         if (k != j) {
            f1 = ffNeg(ffExtract(xk1,j));
            ffAddMulRow(xk1,xj1,f1, noc);
            ffAddMulRow(xk2,xj2,f1, noc);
         }
         ffStepPtr(&xk1, noc);
         ffStepPtr(&xk2, noc);
      }
      ffStepPtr(&xj1, noc);
      ffStepPtr(&xj2, noc);
   }
}


/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Calculates and returns the inverse of a matrix. @p mat must be a non-singular square matrix,
/// otherwise the program is aborted with an error message.
/// The return value is a independent matrix, the original matrix remains unchanged.

Matrix_t *matInverse(const Matrix_t *mat)
{

   matValidate(MTX_HERE, mat);
   if (mat->nor != mat->noc) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_NOTSQUARE);
   }
   const int dim = mat->nor;
   Matrix_t *dest = matId(mat->field, dim);

   // Copy matrix into workspace
   PTR tmp = ffAlloc(dim, dim);
   memcpy(tmp,mat->data,ffSize(dim, dim));

   // Inversion
   zmatinv(tmp,dest->data, dim);
   ffFree(tmp);

   return dest;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
