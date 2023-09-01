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
   uint32_t fileHeader[3];
   sysRead32(f,fileHeader,3);
   if (fileHeader[0] != MTX_TYPE_INTMATRIX)
      mtxAbort(MTX_HERE,"Not an integer matrix (type 0x%lx)", (unsigned long) fileHeader[0]);
   IntMatrix_t *m = imatAlloc(fileHeader[1],fileHeader[2]);
   sysRead32(f,m->Data,m->Nor * m->Noc);
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
