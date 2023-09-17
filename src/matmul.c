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
   if ((src->Field != dest->Field) || (src->Nor != dest->Noc)) {
      mtxAbort(MTX_HERE,"Can't multiply %dx%d/GF(%d) by %dx%d/GF(%d): %s",
                 dest->Nor,dest->Noc,dest->Field,src->Nor,src->Noc,src->Field,
                 MTX_ERR_INCOMPAT);
   }
#endif

   // matrix multiplication
   ffSetField(src->Field);
   result = tmp = ffAlloc(dest->Nor, src->Noc);
   x = dest->Data;
   for (i = 0; i < dest->Nor; ++i) {
      ffMapRow(x,src->Data,src->Nor,src->Noc,tmp);
      ffStepPtr(&tmp, src->Noc);
      ffStepPtr(&x, dest->Noc);
   }
   sysFree(dest->Data);
   dest->Data = result;
   dest->Noc = src->Noc;

   mat_DeletePivotTable(dest);

   return dest;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
