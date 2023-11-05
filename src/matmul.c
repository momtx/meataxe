////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Matrix multiplication
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Multiply matrices
/// This function multiplies @em dest from the right by @em src.
/// The matrices must be compatible for multiplication, i.e. they must be over
/// the same field, and the number of columns of @em dest must be equal to the
/// number of rows of @em src.
/// The result of the multiplication is stored in @em dest, overwriting the
/// original contents.
/// @see matPower()
/// @param dest Left factor and result.
/// @param src Right factor.
/// @return The function returns @em dest.

Matrix_t *matMul(Matrix_t *dest, const Matrix_t *src)
{
   PTR x, tmp, result;
   long i;

   // check arguments
#ifdef MTX_DEBUG
   matValidate(MTX_HERE, src);
   matValidate(MTX_HERE, dest);
   if ((src->field != dest->field) || (src->nor != dest->noc)) {
      mtxAbort(MTX_HERE,"Can't multiply %dx%d/GF(%d) by %dx%d/GF(%d): %s",
                 dest->nor,dest->noc,dest->field,src->nor,src->noc,src->field,
                 MTX_ERR_INCOMPAT);
   }
#endif

   // matrix multiplication
   ffSetField(src->field);
   result = tmp = ffAlloc(dest->nor, src->noc);
   x = dest->data;
   for (i = 0; i < dest->nor; ++i) {
      ffMapRow(tmp, x,src->data,src->nor,src->noc);
      ffStepPtr(&tmp, src->noc);
      ffStepPtr(&x, dest->noc);
   }
   sysFree(dest->data);
   dest->data = result;
   dest->noc = src->noc;

   mat_DeletePivotTable(dest);

   return dest;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
