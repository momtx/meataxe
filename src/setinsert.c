////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Insert an element into a set
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


const int BlockSize = 5;    /* Allocation unit */

/// @addtogroup intset
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Insert an element into a set.
/// @param set Pointer to the set.
/// @param elem Number to insert.
/// @return 0 on success, -1 on error.

int setInsert(Set_t *set, long elem)
{
   register int i, k;
   register long *d;

   setValidate(MTX_HERE, set);

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
         mtxAbort(MTX_HERE,"Cannot grow set");
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
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
