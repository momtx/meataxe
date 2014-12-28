////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Multiplication of permutations
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

/// @addtogroup perm
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Multiply permutations.
/// This function multiplies @em dest from the right by @em src. Both
/// permutations must have the same degree.
/// @param dest Pointer to the first permutation.
/// @param src Pointer to the second permutation.
/// @return @em dest, or 0 on error.

Perm_t *PermMul(Perm_t *dest, const Perm_t *src)
{
   register long i;
   register long *d, *s;

   // Check arguments
   if (!PermIsValid(dest) || !PermIsValid(src)) {
      return NULL;
   }
   if (dest->Degree != src->Degree) {
      MTX_ERROR1("%E",MTX_ERR_INCOMPAT);
      return NULL;
   }

   // Multiply
   d = dest->Data;
   s = src->Data;
   for (i = dest->Degree; i > 0; --i) {
      *d = s[*d];
      ++d;
   }
   return dest;
}


/// @}