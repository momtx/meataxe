////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Read an integer matrix from a file
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup imat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read an integer matrix from a file.
/// This function reads an integer matrix from a file.
/// @see imatLoad()
/// @param f File to read from.
/// @return Pointer to the matrix, or 0 on error.

IntMatrix_t *imatRead(FILE *f)
{
   IntMatrix_t *m;
   long hdr[3];

   if (sysReadLong32(f,hdr,3) != 3) {
      mtxAbort(MTX_HERE,"Cannot read header");
      return 0;
   }
   if (hdr[0] != -8) {  // HACK: T_IMAT in 2.3 
      mtxAbort(MTX_HERE,"Not an integer matrix");
      return 0;
   }
   m = imatAlloc(hdr[1],hdr[2]);
   if (m == 0) {
      return 0;
   }
   if (sysReadLong32(f,m->Data,m->Nor * m->Noc) != m->Nor * m->Noc) {
      imatFree(m);
      return 0;
   }
   return m;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read an integer matrix from a file.
/// This function opens a file, reads a single integer matrix, and closes
/// the file. To read more than one matrix from a file, use imatRead().
/// @param fn File name.
/// @return Pointer to the matrix, or 0 on error.

IntMatrix_t *imatLoad(const char *fn)
{
   FILE *f;
   IntMatrix_t *m;

   if ((f = sysFopen(fn,"rb")) == 0) {
      return 0;
   }
   m = imatRead(f);
   fclose(f);
   return m;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
