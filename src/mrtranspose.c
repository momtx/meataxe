////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Matrix representations, transpose
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

/// @addtogroup mrep
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Transpose a Representation.
/// This function transposes a matrix representation. A new representation
/// is created, and the original representation is not changed.
/// @param rep Matrix representation.
/// @return Pointer to the transposed representation, or 0 on error.

MatRep_t *MrTransposed(const MatRep_t *rep)
{
   Matrix_t **tr;
   MatRep_t *tr_rep;
   int i;

   /* Check arguments
      --------------- */
   if (!MrIsValid(rep)) {
      MTX_ERROR1("rep: %E",MTX_ERR_BADARG);
      return NULL;
   }

   /* Transpose the generators
      ------------------------ */
   tr = NALLOC(Matrix_t *,rep->NGen);
   if (tr == NULL) {
      MTX_ERROR("Cannot allocate buffer");
      return NULL;
   }
   for (i = 0; i < rep->NGen; ++i) {
      tr[i] = MatTransposed(rep->Gen[i]);
      if (tr[i] == NULL) {
         while (--i > 0) {
            MatFree(tr[i]);
         }
         SysFree(tr);
         MTX_ERROR("Cannot transpose generator");
         return NULL;
      }
   }

   /* Make the new representation
      --------------------------- */
   tr_rep = MrAlloc(rep->NGen,tr,0);
   if (tr_rep == NULL) {
      for (i = 0; i < rep->NGen; ++i) {
         MatFree(tr[i]);
      }
      SysFree(tr);
      return NULL;
   }

   SysFree(tr);
   return tr_rep;
}


/// @}
