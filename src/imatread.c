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

/// Reads an integer matrix from a file. See also @ref imatLoad.

IntMatrix_t* imatRead(MtxFile_t* file)
{
   mfReadHeader(file);
   if (mfObjectType(file) != MTX_TYPE_INTMATRIX) {
      mtxAbort(MTX_HERE, "%s: unexpected object type 0x%lx (expected integer matrix)",
               file->name, (unsigned long) file->header[0]);
   }
   IntMatrix_t* m = imatAlloc(file->header[1], file->header[2]);
   mfRead32(file, m->data, m->nor * m->noc);
   return m;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Read an integer matrix from a file.
/// This function opens a file, reads a single integer matrix, and closes the file.
/// See also use @ref imatRead.
/// @param fn File name.
/// @return Pointer to the matrix.

IntMatrix_t *imatLoad(const char *fn)
{
   MtxFile_t* file = mfOpen(fn);
   IntMatrix_t* m = imatRead(file);
   mfClose(file);
   return m;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
