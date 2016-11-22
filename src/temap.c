////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Map under tensor product.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
//   Local data

MTX_DEFINE_FILE_INFO

/// @addtogroup tp
/// @{

/// Map under tensor product.
/// This function applies the tensor product of two matrices to one or more
/// vectors. The same calculation could be done with MatMul() and
/// MatTensor(), but this function is usually faster and uses less memory,
/// because it does not calculate the full tensor product of a⊗b.
/// @see VectorToMatrix() MatrixToVector()
/// @param vec Vectors to map.
/// @param a Left matrix.
/// @param b Right matrix.
/// @return Image of @a vec under @a a⊗@a b, or 0 on error.

Matrix_t *TensorMap(Matrix_t *vec, const Matrix_t *a, const Matrix_t *b)
{

   // Check the arguments
   if (!MatIsValid(vec)) {
      MTX_ERROR1("vec: %E",MTX_ERR_BADARG);
      return NULL;
   }
   if (!MatIsValid(a)) {
      MTX_ERROR1("a: %E",MTX_ERR_BADARG);
      return NULL;
   }
   if (!MatIsValid(b)) {
      MTX_ERROR1("b: %E",MTX_ERR_BADARG);
      return NULL;
   }
   if ((a->Field != b->Field) || (b->Field != vec->Field) ||
       (vec->Noc != a->Nor * b->Nor)) {
      MTX_ERROR1("%E",MTX_ERR_INCOMPAT);
      return NULL;
   }

   // Calculate the result
   Matrix_t *result = MatAlloc(vec->Field,vec->Nor,a->Noc * b->Noc);
   if (result == NULL) {
      return NULL;
   }
   for (int i = 0; i < vec->Nor; ++i) {
      Matrix_t *tmp = MatTransposed(a);
      if (tmp == NULL) {
         MatFree(result);
         return NULL;
      }

      Matrix_t *v = VectorToMatrix(vec,i,b->Nor);
      if (v == NULL) {
         MTX_ERROR("Conversion failed");
         MatFree(result);
         return NULL;
      }
      if (MatMul(tmp,v) == NULL) {
         MatFree(result);
         return NULL;
      }
      MatFree(v);
      if (MatMul(tmp,b) == NULL || MatrixToVector(tmp,result,i) != 0) {
         MatFree(result);
         MTX_ERROR("Conversion failed");
         return NULL;
      }
      MatFree(tmp);
   }
   return result;
}


/// @}
