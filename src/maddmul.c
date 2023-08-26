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
   if ((dest->Field != src->Field) || (dest->Nor != src->Nor) ||
       (dest->Noc != src->Noc)) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
      return NULL;
   }

   /* Handle special cases
      -------------------- */
   if (coeff == FF_ONE) {
      matAdd(dest,src);
   } else if (coeff == FF_ZERO) {
   } else {
      /* Add multiple
         ------------ */
      PTR dp = dest->Data, sp = src->Data;
      int n;
      ffSetField(src->Field);
      ffSetNoc(src->Noc);
      for (n = src->Nor; n > 0; --n) {
         ffAddMulRow(dp,sp,coeff);
         ffStepPtr(&dp, src->Noc);
         ffStepPtr(&sp, src->Noc);
      }
   }
   return dest;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
