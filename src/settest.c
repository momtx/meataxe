////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check if a set contains an element
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

/// @addtogroup intset
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Check if a number is in a set.
/// @param set Pointer to the set.
/// @param elem Number to check.
/// @return 1 if the set contains @a elem, 0 if not.

int SetContains(const Set_t *set, long elem)
{
   int i;
   long *l;

   if (!SetIsValid(set)) {
      MTX_ERROR1("%E",MTX_ERR_BADARG);
      return -1;
   }
   l = set->Data;
   for (i = set->Size; i > 0 && *l < elem; --i, ++l) {
   }
   return i > 0 && *l == elem;
}


/// @}
