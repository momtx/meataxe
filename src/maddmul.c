////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Add multiple of a matrix
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Add a multiple of a matrix.
/// This function adds a multiple of a matrix to another matrix. The
/// matrices must be compatible for addition.
/// |matAddMul()| handles special cases (|coeff| equals 0 or 1) in an
/// intelligent way, so there is no need for the caller to do this.
/// @param dest Matrix to add to.
/// @param src Matrix to add.
/// @param coeff Coefficient.
/// @return The function returns |dest|, or |NULL| on error.

Matrix_t *matAddMul(Matrix_t *dest, const Matrix_t *src, FEL coeff)
{
   // Check arguments
   matValidate(MTX_HERE,src);
   matValidate(MTX_HERE,dest);
   if ((dest->field != src->field) || (dest->nor != src->nor) ||
       (dest->noc != src->noc)) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
      return NULL;
   }

   if (coeff == FF_ONE) {
      matAdd(dest,src);
   }
   else if (coeff != FF_ZERO) {
      PTR dp = dest->data, sp = src->data;
      int n;
      ffSetField(src->field);
      for (n = src->nor; n > 0; --n) {
         ffAddMulRow(dp,sp,coeff,src->noc);
         ffStepPtr(&dp, src->noc);
         ffStepPtr(&sp, src->noc);
      }
   }
   return dest;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
