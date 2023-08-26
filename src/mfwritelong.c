////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Write long integers
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup mf
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write long integers to a file.
/// This function writes |count| long integers from buffer into a data file.
/// If necessary, the data is converted from machine independent format into
/// the format needed by the platform. See |sysWriteLong()| for details.
/// @param f Pointer to the file.
/// @param buf Data buffer.
/// @param count Number of integers to write.
/// @return Number of integers that were actually written. Any value other than count
///    indicates an error.

int mfWriteLong(MtxFile_t *f, const long *buf, int count)
{
   int rc;
   if (!mfIsValid(f)) {
      return -1;
   }
   rc = sysWriteLong32(f->File,buf,count);
   if (rc != count) {
      mtxAbort(MTX_HERE,"%s: write failed",f->Name);
   }
   return rc;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
