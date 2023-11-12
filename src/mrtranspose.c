////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Matrix representations, transpose
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup mrep
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Transpose a Representation.
/// This function transposes a matrix representation. A new representation
/// is created, and the original representation is not changed.
/// @param rep Matrix representation.
/// @return Pointer to the transposed representation, or 0 on error.

MatRep_t *mrTransposed(const MatRep_t *rep)
{
   mrValidate(MTX_HERE,rep);

   Matrix_t **tr = NALLOC(Matrix_t *,rep->NGen);
   for (int i = 0; i < rep->NGen; ++i) {
      tr[i] = matTransposed(rep->Gen[i]);
   }

   MatRep_t *tr_rep = mrAlloc(rep->NGen,tr,0);

   sysFree(tr);
   return tr_rep;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
