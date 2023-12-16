////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Find the (kernel-)stable power of a matrix.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
//   Local data


/// @addtogroup algo
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Stable power of a matrix.
/// This function takes a square matrix M and finds an integer n>0 such that
/// ker(M<sup>n</sup>) = ker(M<sup>n+1</sup>).
/// @a ker must be a pointer to a variable of type Matrix_t*,
/// where the stable kernel will be stored. Both @a pwr and @a ker may be
/// NULL if the corresponding information is not needed.
///
/// Note that the number $n$ found by StablePower_() is not guararanteed
/// to be minimal. In fact, n will always be a power of two since the
/// function only examines matrices of the form M<sup>2<sup>k</sup></sup>.
///
/// This function modifies the matrix. To avoid this, use StablePower().
/// @param mat The matrix.
/// @param pwr Stable power.
/// @param ker Kernel of the stable power.
/// @return 0 on success, -1 on error.

int StablePower_(Matrix_t *mat, int *pwr, Matrix_t **ker)
{
   // check the arguments.
   matValidate(MTX_HERE, mat);
   if (mat->nor != mat->noc) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_NOTSQUARE);
   }

   // calculate the stable power
   int p = 1;
   Matrix_t *k1 = matNullSpace(mat);
   if (k1 == NULL) {
       return -1;
   }
   if (matMul(mat,mat) == NULL) {
       matFree(k1);
       return -1;
   }
   Matrix_t *k2 = matNullSpace(mat);
   if (k2 == NULL) {
       matFree(k1);
       return -1;
   }

   while (k2->nor > k1->nor) {
      p *= 2;
      matFree(k1);
      k1 = k2;
      if (matMul(mat,mat) == NULL) {
         matFree(k1);
         return -1;
      }
      k2 = matNullSpace(mat);
      if (k2 == NULL) {
         matFree(k1);
         return -1;
      }
   }
   matFree(k2);

   // return the result
   if (ker != NULL) {
      *ker = k1;
   } else {
      matFree(k1);
   }
   if (pwr != NULL) {
      *pwr = p;
   }

   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Stable power of a matrix.
/// This function works like |StablePower_()|, but it does not modify
/// the matrix  in |mat|. This means, of course, that a temporary copy
/// of the matrix is created.
/// @param mat The matrix.
/// @param pwr Stable power.
/// @param ker Kernel of the stable power.
/// @return $0$ on success, $-1$ on error.

int StablePower(const Matrix_t *mat, int *pwr, Matrix_t **ker)
{
   int rc;
   Matrix_t *tmp;

   tmp = matDup(mat);
   if (tmp == NULL) {
      mtxAbort(MTX_HERE,"mat: %s",MTX_ERR_BADARG);
      return -1;
   }
   rc = StablePower_(tmp,pwr,ker);
   matFree(tmp);
   return rc;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
