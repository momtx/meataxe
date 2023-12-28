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

Matrix_t* matInsert_(Matrix_t* mat, const Poly_t* pol)
{
   int i;
   int l;
   PTR v;
   FEL f;

   matValidate(MTX_HERE, mat);
   polValidate(MTX_HERE, pol);
   const uint32_t nor = mat->nor;
   if (nor != mat->noc || mat->field != pol->field) {
      mtxAbort(MTX_HERE, "%s", MTX_ERR_NOTSQUARE);
   }

   ffSetField(mat->field);

   // Special case: p(x) = 0
   if (pol->degree == -1) {
      for (l = 0, v = mat->data; l < nor; ffStepPtr(&v, nor), ++l) {
         ffMulRow(v, FF_ZERO, nor);
      }
      return mat;
   }

   // Special case: deg(p) = 0
   if (pol->degree == 0) {
      for (l = 0, v = mat->data; l < nor; ffStepPtr(&v, nor), ++l) {
         ffMulRow(v, FF_ZERO, nor);
         ffInsert(v, l, pol->data[0]);
      }
      return mat;
   }

   // Evaluate p(A)
   Matrix_t* x = NULL;
   if (pol->degree > 1) {
      x = matDup(mat);
   }
   if ((f = pol->data[pol->degree]) != FF_ONE) {
      for (l = nor, v = mat->data; l > 0; --l, ffStepPtr(&v, nor)) {
         ffMulRow(v, f, nor);
      }
   }
   for (i = pol->degree - 1; i >= 0; --i) {
      if ((f = pol->data[i]) != FF_ZERO) {
         for (l = 0, v = mat->data; l < nor; ++l, ffStepPtr(&v, nor)) {
            ffInsert(v, l, ffAdd(ffExtract(v, l), f));
         }
      }
      if (i > 0) {
         matMul(mat, x);
      }
   }
   if (x != NULL) {
      matFree(x);
   }
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
   uint32_t nor;
   int l;
   PTR v;
   FEL f;

   matValidate(MTX_HERE, mat);
   polValidate(MTX_HERE, pol);
   if ((nor = mat->nor) != mat->noc) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_NOTSQUARE);
      return NULL;
   }
   if (mat->field != pol->field) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
      return NULL;
   }

   // Special case: p = 0
   if (pol->degree == -1) {
      return matAlloc(mat->field,nor,nor);
   }

   // Special case: deg(p) = 0
   if (pol->degree == 0) {
      x = matAlloc(mat->field,nor,nor);
      if (x == NULL) {
	  return NULL;
      }
      for (l = 0, v = x->data; l < nor; ++l, ffStepPtr(&v, nor)) {
         ffInsert(v,l,pol->data[0]);
      }
      return x;
   }

   // Evaluate p(A)
   x = matDup(mat);
   if (x == NULL) {
      return NULL;
   }
   if ((f = pol->data[pol->degree]) != FF_ONE) {
      for (l = nor, v = x->data; l > 0; --l, ffStepPtr(&v, nor)) {
         ffMulRow(v,f, nor);
      }
   }
   for (i = pol->degree - 1; i >= 0; --i) {
      if ((f = pol->data[i]) != FF_ZERO) {
         for (l = 0, v = x->data; l < nor; ++l, ffStepPtr(&v, nor)) {
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
