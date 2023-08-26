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
   Matrix_t **tr;
   MatRep_t *tr_rep;
   int i;

   /* Check arguments
      --------------- */
   if (!mrIsValid(rep)) {
      mtxAbort(MTX_HERE,"rep: %s",MTX_ERR_BADARG);
      return NULL;
   }

   /* Transpose the generators
      ------------------------ */
   tr = NALLOC(Matrix_t *,rep->NGen);
   if (tr == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate buffer");
      return NULL;
   }
   for (i = 0; i < rep->NGen; ++i) {
      tr[i] = matTransposed(rep->Gen[i]);
      if (tr[i] == NULL) {
         while (--i > 0) {
            matFree(tr[i]);
         }
         sysFree(tr);
         mtxAbort(MTX_HERE,"Cannot transpose generator");
         return NULL;
      }
   }

   /* Make the new representation
      --------------------------- */
   tr_rep = mrAlloc(rep->NGen,tr,0);
   if (tr_rep == NULL) {
      for (i = 0; i < rep->NGen; ++i) {
         matFree(tr[i]);
      }
      sysFree(tr);
      return NULL;
   }

   sysFree(tr);
   return tr_rep;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
