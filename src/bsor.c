////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Bit Strings, "or" operation
//
// (C) Copyright 1998-2014 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

/// @addtogroup bs
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Union of two bit strings.
/// This function computes the union of two bit strings, i.e., the bitwise logical "or".
/// The result is stored in @em dest and overwrites the previous value. Both bit
/// strings must have the same size.
/// @return 0 on success, -1 on error.

int BsOr(BitString_t *dest, const BitString_t *src)
{
   register int i;
   register unsigned long *dp;
   register const unsigned long *sp;

   // check arguments
   if (!BsIsValid(dest)) {
      MTX_ERROR1("dest: %E",MTX_ERR_BADARG);
      return -1;
   }
   if (!BsIsValid(src)) {
      MTX_ERROR1("src: %E",MTX_ERR_BADARG);
      return -1;
   }
   if (dest->Size != src->Size) {
      MTX_ERROR1("%E",MTX_ERR_INCOMPAT);
      return -1;
   }

   // OR operation
   dp = (unsigned long *) dest->Data;
   sp = (unsigned long const *) src->Data;
   for (i = src->BufSize; i > 0; --i) {
      *dp++ |= *sp++;
   }

   return 0;
}


/// @}
