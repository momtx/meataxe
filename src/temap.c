////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Map under tensor product.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
//   Local data


/// @addtogroup tp
/// @{

/// Map under tensor product.
/// This function applies the tensor product of two matrices to one or more
/// vectors. The same calculation could be done with matMul() and
/// matTensor(), but this function is usually faster and uses less memory,
/// because it does not calculate the full tensor product of a⊗b.
/// @see VectorToMatrix() MatrixToVector()
/// @param vec Vectors to map.
/// @param a Left matrix.
/// @param b Right matrix.
/// @return Image of @p vec under @p a⊗@p b, or 0 on error.

Matrix_t *TensorMap(Matrix_t *vec, const Matrix_t *a, const Matrix_t *b)
{

   // Check the arguments
   matValidate(MTX_HERE, vec);
   matValidate(MTX_HERE, a);
   matValidate(MTX_HERE, b);
   if ((a->field != b->field) || (b->field != vec->field) ||
       (vec->noc != a->nor * b->nor)) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
      return NULL;
   }

   // Calculate the result
   Matrix_t *result = matAlloc(vec->field,vec->nor,a->noc * b->noc);
   if (result == NULL) {
      return NULL;
   }
   for (int i = 0; i < vec->nor; ++i) {
      Matrix_t *tmp = matTransposed(a);
      if (tmp == NULL) {
         matFree(result);
         return NULL;
      }

      Matrix_t *v = VectorToMatrix(vec,i,b->nor);
      if (v == NULL) {
         mtxAbort(MTX_HERE,"Conversion failed");
         matFree(result);
         return NULL;
      }
      if (matMul(tmp,v) == NULL) {
         matFree(result);
         return NULL;
      }
      matFree(v);
      if (matMul(tmp,b) == NULL || MatrixToVector(tmp,result,i) != 0) {
         matFree(result);
         mtxAbort(MTX_HERE,"Conversion failed");
         return NULL;
      }
      matFree(tmp);
   }
   return result;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
