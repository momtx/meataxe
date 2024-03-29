////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Power of a matrix
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Calculate the n-th power of a matrix using the binary method.
/// @param n The exponent.
/// @param inp Input matrix (dim x dim)
/// @param out Output matrix (dim x dim)
/// @param tmp2 Workspace (1xdim)
/// @param dim Matrix size

static void matpwr_(long n, PTR inp, PTR out, PTR tmp2, int dim)
{
   PTR x, y;
   long i;
   int first = 1;

   while (n > 0) {
      if (n % 2 == 1) {
         if (first) {
            memcpy(out,inp,ffSize(dim, dim));
            first = 0;
         } else {
            x = out;
            for (i = 0; i < dim; ++i) {
               ffMapRow(tmp2, x,inp,dim,dim);
               ffCopyRow(x,tmp2, dim);
               ffStepPtr(&x, dim);
            }
         }
      }
      if (n == 1) {
         break;
      }
      x = inp;
      y = tmp2;
      for (i = 0; i < dim; ++i) {
         ffMapRow(y, x,inp,dim, dim);
         ffStepPtr(&x, dim);
         ffStepPtr(&y, dim);
      }
      memcpy(inp,tmp2,ffSize(dim, dim));
      n /= 2;
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Power of a matrix.
/// This function calculates the n-th power of a matrix, using the binary
/// method. This method is generally faster than multiplying the matrix $n$
/// times by itself. On the other hand, a third matrix is temporarily created
/// in addition to the original matrix and the result matrix.
/// The cases n=0 and n=1 are handled separately, avoiding unnecessary
/// memory allocation and calculation.
///
/// Negative exponents are not allowed. To calculate a negative power, you
/// must first invert the matrix with matInverse() and then call matPower()
/// with the inverted matrix and a positive exponent.
/// @param mat Pointer to the matrix.
/// @param n Exponent.
/// @return n-th power of mat, or NULL on error.

Matrix_t *matPower(const Matrix_t *mat, long n)
{
   Matrix_t *result;
   PTR tmp, tmp2;

   /* Check the arguments
      ------------------- */
   matValidate(MTX_HERE, mat);
   if (mat->nor != mat->noc) {
      mtxAbort(MTX_HERE,"matPower(): %s",MTX_ERR_NOTSQUARE);
      return NULL;
   }

   /* Handle special cases n = 0 and n = 1
      ------------------------------------ */
   if (n == 0) {
      return matId(mat->field,mat->nor);
   } else if (n == 1) {
      return matDup(mat);
   }

   ffSetField(mat->field);
   tmp = ffAlloc(mat->noc, mat->noc);
   if (tmp == NULL) {
       return NULL;
   }
   memcpy(tmp,mat->data,ffSize(mat->noc, mat->noc));
   tmp2 = ffAlloc(mat->noc, mat->noc);
   if (tmp2 == NULL) {
       sysFree(tmp);
       return NULL;
   }
   result = matAlloc(mat->field,mat->nor,mat->noc);
   if (result != NULL) {
      matpwr_(n, tmp, result->data, tmp2, mat->noc);
   }
   sysFree(tmp);
   sysFree(tmp2);
   return result;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
