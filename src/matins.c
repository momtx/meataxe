////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Insert a matrix into a polynomial
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Insert a matrix into a polynomial
/// Given a square matrix A and a polynomial p over the same field, this functions
/// calculates p(A). Unlike matInsert() this function is destructive. The result
/// is stored in the original matrix and the old value is lost.
/// @param mat Pointer to the matrix.
/// @param pol Pointer to the polynomial.
/// @return The function returns @em mat, or 0 on error.

Matrix_t *matInsert_(Matrix_t *mat, const Poly_t *pol)
{
   Matrix_t *x = NULL;
   int i;
   int nor;
   int l;
   PTR v;
   FEL f;

   // Check the arguments
   matValidate(MTX_HERE, mat);
   polValidate(MTX_HERE, pol);
   if ((nor = mat->Nor) != mat->Noc) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_NOTSQUARE);
      return NULL;
   }
   if (mat->Field != pol->Field) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
      return NULL;
   }

   ffSetField(mat->Field);

   // Special case: p(x) = 0
   if (pol->Degree == -1) {
      for (l = 0, v = mat->Data; l < nor; ffStepPtr(&v, nor), ++l) {
         ffMulRow(v,FF_ZERO, nor);
      }
      return mat;
   }

   // Special case: deg(p) = 0
   if (pol->Degree == 0) {
      for (l = 0, v = mat->Data; l < nor; ffStepPtr(&v, nor), ++l) {
         ffMulRow(v,FF_ZERO, nor);
         ffInsert(v,l,pol->Data[0]);
      }
      return mat;
   }

   // Evaluate p(A)
   if (pol->Degree > 1) {
      x = matDup(mat);
      if (x == NULL) {
	  return NULL;
      }
   }
   if ((f = pol->Data[pol->Degree]) != FF_ONE) {
      for (l = nor, v = mat->Data; l > 0; --l, ffStepPtr(&v, nor)) {
         ffMulRow(v,f, nor);
      }
   }
   for (i = pol->Degree - 1; i >= 0; --i) {
      if ((f = pol->Data[i]) != FF_ZERO) {
         for (l = 0, v = mat->Data; l < nor; ++l, ffStepPtr(&v, nor)) {
            ffInsert(v,l,ffAdd(ffExtract(v,l),f));
         }
      }
      if (i > 0) {
         matMul(mat,x);
      }
   }
   if (pol->Degree > 1) { matFree(x); }
   return mat;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Insert a matrix into a polynomial
/// Given a square matrix A and a polynomial p over the same field, this functions
/// calculates p(A). Unlike matInsert_() this function returns a new matrix and does
/// not modify the original matrix.
/// @param mat Pointer to the matrix.
/// @param pol Pointer to the polynomial.
/// @return @em pol(@em mat), or 0 on error.

Matrix_t *matInsert(const Matrix_t *mat, const Poly_t *pol)
{
   Matrix_t *x;
   int i;
   int nor;
   int l;
   PTR v;
   FEL f;

   matValidate(MTX_HERE, mat);
   polValidate(MTX_HERE, pol);
   if ((nor = mat->Nor) != mat->Noc) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_NOTSQUARE);
      return NULL;
   }
   if (mat->Field != pol->Field) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
      return NULL;
   }

   // Special case: p = 0
   if (pol->Degree == -1) {
      return matAlloc(mat->Field,nor,nor);
   }

   // Special case: deg(p) = 0
   if (pol->Degree == 0) {
      x = matAlloc(mat->Field,nor,nor);
      if (x == NULL) {
	  return NULL;
      }
      for (l = 0, v = x->Data; l < nor; ++l, ffStepPtr(&v, nor)) {
         ffInsert(v,l,pol->Data[0]);
      }
      return x;
   }

   // Evaluate p(A)
   x = matDup(mat);
   if (x == NULL) {
      return NULL;
   }
   if ((f = pol->Data[pol->Degree]) != FF_ONE) {
      for (l = nor, v = x->Data; l > 0; --l, ffStepPtr(&v, nor)) {
         ffMulRow(v,f, nor);
      }
   }
   for (i = pol->Degree - 1; i >= 0; --i) {
      if ((f = pol->Data[i]) != FF_ZERO) {
         for (l = 0, v = x->Data; l < nor; ++l, ffStepPtr(&v, nor)) {
            ffInsert(v,l,ffAdd(ffExtract(v,l),f));
         }
      }
      if (i > 0) {
         matMul(x,mat);
      }
   }
   return x;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
