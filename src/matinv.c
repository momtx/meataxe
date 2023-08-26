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

   MTX_ASSERT(ffNoc == noc, );
   // initialize result with identity matrix
   for (j = 0, xj1 = result; j < noc; ++j, ffStepPtr(&xj1, noc)) {
      ffMulRow(xj1,FF_ZERO);
      ffInsert(xj1,j,FF_ONE);
   }

   // matrix inversion
   xj1 = mat;
   xj2 = result;
   for (j = 0; j < noc; ++j) {

      for (xk1 = xj1, k = j;
           k < ffNoc && (f1 = ffExtract(xk1,j)) == FF_ZERO;
           ++k, ffStepPtr(&xk1, noc)) {
      }
      if (f1 == FF_ZERO) {
         mtxAbort(MTX_HERE,"%s",MTX_ERR_DIV0);
      }
      if (k > j) {      /* Swap rows */
         ffSwapRows(xk1,xj1);
/*	    xk2 = ffGetPtr(xj2,k-j,noc);*/
         xk2 = (PTR)((char *)xj2 + (k - j) * ffRowSize(noc));
         ffSwapRows(xk2,xj2);
      }
      f2 = ffInv(f1);
      ffMulRow(xj1,f2);
      ffMulRow(xj2,f2);
      xk1 = mat;
      xk2 = result;
      for (k = 0; k < ffNoc; ++k) {
         if (k != j) {
            f1 = ffNeg(ffExtract(xk1,j));
            ffAddMulRow(xk1,xj1,f1);
            ffAddMulRow(xk2,xj2,f1);
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
/// Calculates and returns the inverse of a matrix. @a mat must be a non-singular square matrix,
/// otherwise the program is aborted with an error message.
/// The return value is a independent matrix, the original matrix remains unchanged.

Matrix_t *matInverse(const Matrix_t *mat)
{
   PTR tmp = NULL;      // workspace

   matValidate(MTX_HERE, mat);
   if (mat->Nor != mat->Noc) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_NOTSQUARE);
   }
   const int dim = mat->Nor;
   Matrix_t *dest = matId(mat->Field, dim);

   // Copy matrix into workspace
   tmp = ffAlloc(dim, dim);
   memcpy(tmp,mat->Data,ffSize(dim, dim));

   // Inversion
   zmatinv(tmp,dest->Data, dim);
   return dest;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
