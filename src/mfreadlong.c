////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Read long integers
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup mf
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read long integers from a file.
/// This function reads |count| long integers from a data file into a buffer.
/// If necessary, the data is converted from machine independent format into
/// the format needed by the platform. See |sysReadLong()| for details.
/// @param f Pointer to the file.
/// @param buf Data buffer.
/// @param count Number of integers to read.
/// @return Number of integers that were actually read.
///    Any value other than |count| indicates an error.

int mfReadLong(MtxFile_t *f, long *buf, int count)
{
   int rc;
   if (!mfIsValid(f)) {
      return -1;
   }
   rc = sysReadLong32(f->File,buf,count);
   if (rc < 0) {
      mtxAbort(MTX_HERE,"%s: read failed",f->Name);
   }
   return rc;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
