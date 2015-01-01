////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Insert an element into a set
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

const int BlockSize = 5;    /* Allocation unit */

/// @addtogroup intset
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Insert an element into a set.
/// @param set Pointer to the set.
/// @param elem Number to insert.
/// @return 0 on success, -1 on error.

int SetInsert(Set_t *set, long elem)
{
   register int i, k;
   register long *d;

   if (!SetIsValid(set)) {
      MTX_ERROR1("set: %E",MTX_ERR_BADARG);
      return -1;
   }

   // Find the position of <elem> in the data buffer.
   // If <elem> is already containd in the set, exit.
   // TODO: binary search
   d = set->Data;
   for (i = 0; i < set->Size && d[i] < elem; ++i) {
   }
   if ((i < set->Size) && (d[i] == elem)) {
      return 0;
   }

   // extend the buffer if necessary
   if (set->Size >= set->BufSize) {
      int newmax = set->BufSize + BlockSize;
      long *newbuf = NREALLOC(set->Data,long,newmax);
      if (newbuf == NULL) {
         MTX_ERROR("Cannot grow set");
         return -1;
      }
      set->BufSize = newmax;
      d = set->Data = newbuf;
   }

   // insert new element
   for (k = set->Size - 1; k >= i; --k) {
      d[k + 1] = d[k];
   }
   d[i] = elem;
   ++set->Size;

   return 0;
}


/// @}
