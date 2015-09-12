////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Duplicate a permutation
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

/// @addtogroup perm
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Duplicate a permutation.
/// This function creates a copy of an existing permutation.
/// @param src Pointer to the permutation.
/// @return Pointer to a copy of @em src or 0 on error.

Perm_t *PermDup(const Perm_t *src)
{
   Perm_t *p;

   if (!PermIsValid(src)) {
      MTX_ERROR1("src: %E",MTX_ERR_BADARG);
      return NULL;
   }
   p = PermAlloc(src->Degree);
   if (p == NULL) {
      MTX_ERROR("Cannot allocate result");
      return NULL;
   }
   memcpy(p->Data,src->Data,(size_t) src->Degree * sizeof(long));
   return p;
}


/// @}
