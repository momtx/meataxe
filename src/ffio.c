////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - I/o for vectors and matrices.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#include <errno.h>
#include <stdlib.h>
#include <string.h>

/// @addtogroup ff
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read Rows
/// This function reads 1 or more rows from a file.
/// @param f Pointer to File.
/// @param buf Pointer to data buffer.
/// @param nor Number of rows to read.
/// @param noc Row size (number of columns).

void ffReadRows(FILE *f, PTR buf, int nor, int noc)
{
   if (noc == 0 || nor == 0) {
      return;
   }
   const size_t rowSizeUsed = ffRowSizeUsed(noc);
   const size_t rowSize = ffRowSize(noc);

   // read rows one by one
   char *b = (char *) buf;
   for (int i = 0; i < nor; ++i) {
      if (fread(b,rowSizeUsed,1,f) != 1) {
         mtxAbort(MTX_HERE,"Error reading row(s) from file: %s", strerror(errno));
      }
      b += rowSize;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write rows
/// This function writes 1 or more rows to a file.
/// The row size must have been set before.
/// @param f Pointer to File.
/// @param buf Pointer to data buffer.
/// @param nor Number of rows to write.
/// @param noc Row size (number of columns).

void ffWriteRows(FILE *f, PTR buf, int nor, int noc)
{
   if (noc == 0 || nor == 0) {
      return;
   }
   size_t rowSizeUsed = ffRowSizeUsed(noc);
   size_t rowSize = ffRowSize(noc);


   const char *b = (const char *) buf;
   for (int i = 0; i < nor; ++i) {
      if (fwrite(b,rowSizeUsed,1,f) != 1) {
         mtxAbort(MTX_HERE,"Error writing row(s) to file: %s", strerror(errno));
      }
      b += rowSize;
   }
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
