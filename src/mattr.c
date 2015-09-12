////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Transpose a matrix
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

/// @addtogroup mat
/// @{

/// Transpose a matrix.
/// @param src Pointer to the matrix.
/// @return Pointer to the transposed matrix or 0 on error.

Matrix_t *MatTransposed(const Matrix_t *src)
{
   PTR s, d;
   long i;
   Matrix_t *dest;

#ifdef DEBUG
   if (!MatIsValid(src)) {
      MTX_ERROR1("src: %E",MTX_ERR_BADARG);
      return NULL;
   }
#endif
   dest = MatAlloc(src->Field,src->Noc,src->Nor);
   if (dest == NULL) {
      MTX_ERROR("Cannot allocate result");
      return NULL;
   }
   d = dest->Data;
   for (i = 0; i < src->Noc; ++i) {
      int k;
      s = src->Data;
      for (k = 0; k < src->Nor; ++k) {
#if defined(DEBUG) && defined(PARANOID)
         FEL f;
         FfSetNoc(src->Noc);
         f = FfExtract(s,i);
         FfSetNoc(src->Nor);
         FfInsert(d,k,f);
#else
         FfInsert(d,k,FfExtract(s,i));
#endif
         s = (PTR)((char*) s + src->RowSize);
      }
      /*d = FfGetPtr(d,1,dest->Noc);*/
      d = (PTR)((char*) d + dest->RowSize);

   }
   return dest;
}


/// @}
