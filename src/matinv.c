////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Matrix inversion
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

////////////////////////////////////////////////////////////////////////////////////////////////////

static int zmatinv(PTR mat, PTR result)
{
   PTR xj1, xj2, xk1, xk2;
   FEL f1 = FF_ZERO, f2;
   long j, k;

   // initialize result with identity matrix
   for (j = 0, xj1 = result; j < FfNoc; ++j, FfStepPtr(&xj1)) {
      FfMulRow(xj1,FF_ZERO);
      FfInsert(xj1,j,FF_ONE);
   }

   // matrix inversion
   xj1 = mat;
   xj2 = result;
   for (j = 0; j < FfNoc; ++j) {

      for (xk1 = xj1, k = j;
           k < FfNoc && (f1 = FfExtract(xk1,j)) == FF_ZERO;
           ++k, FfStepPtr(&xk1)) {
      }
      if (f1 == FF_ZERO) {
         MTX_ERROR1("%E",MTX_ERR_DIV0);
         return -1;
      }
      if (k > j) {      /* Swap rows */
         FfSwapRows(xk1,xj1);
/*	    xk2 = FfGetPtr(xj2,k-j,FfNoc);*/
         xk2 = (PTR)((char *)xj2 + (k - j) * FfCurrentRowSize);
         FfSwapRows(xk2,xj2);
      }
      f2 = FfInv(f1);
      FfMulRow(xj1,f2);
      FfMulRow(xj2,f2);
      xk1 = mat;
      xk2 = result;
      for (k = 0; k < FfNoc; ++k) {
         if (k != j) {
            f1 = FfNeg(FfExtract(xk1,j));
            FfAddMulRow(xk1,xj1,f1);
            FfAddMulRow(xk2,xj2,f1);
         }
         FfStepPtr(&xk1);
         FfStepPtr(&xk2);
      }
      FfStepPtr(&xj1);
      FfStepPtr(&xj2);
   }
   return 0;
}


/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Inverse of a matrix
/// This function calculates the inverse of a matrix. @em mat must be a
/// non-singular square matrix. The inverse matrix is returned in a newly
/// allocated Matrix_t structure, and the original matrix remains unchanged.
/// @param mat Pointer to the matrix.
/// @return Inverse matrix or 0 on error.

Matrix_t *MatInverse(const Matrix_t *mat)
{
   PTR tmp = NULL;      // workspace
   Matrix_t *dest;

   if (!MatIsValid(mat)) {
      return NULL;
   }
   if (mat->Nor != mat->Noc) {
      MTX_ERROR1("%E",MTX_ERR_NOTSQUARE);
      return NULL;
   }
   dest = MatId(mat->Field,mat->Nor);
   if (dest == NULL) {
      return NULL;
   }

   // Copy matrix into workspace
   tmp = FfAlloc(mat->Nor);
   memcpy(tmp,mat->Data,FfCurrentRowSize * mat->Nor);

   // Inversion
   if (zmatinv(tmp,dest->Data) != 0) {
      MatFree(dest);
      return NULL;
   }
   return dest;
}


/// @}