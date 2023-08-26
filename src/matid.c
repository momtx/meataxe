////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Identity matrix
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Identity matrix
/// This function creates an identity matrix with @em nor nows over GF(@em fl).
/// @param fl Field order.
/// @param nor Number of rows.
/// @return Pointer to the matrix, or 0 on error.
/// @see matAlloc

Matrix_t *matId(int fl, int nor)
{
   Matrix_t *m;
   PTR x;
   long i;

   // check arguments
   if ((fl < 2) || (nor < 0)) {
      mtxAbort(MTX_HERE,"Matid(%d,%d): %s",fl,nor,MTX_ERR_BADARG);
      return NULL;
   }

   // Allocate an empty matrix
   m = matAlloc(fl,nor,nor);
   if (m == NULL) {
      return NULL;
   }

   // Set diagonal elements to 1
   for (i = 0, x = m->Data; i < nor; ++i, ffStepPtr(&x, nor)) {
      ffInsert(x,i,FF_ONE);
   }

   return m;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
