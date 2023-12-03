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

/// Writes matrix rows to a binary file.
/// The field must have been set before calling this function, see @ref ffSetField.
/// The function fails and aborts the program if the data could not be written.
///
/// @param file File to write to. The file mut be writable.
/// @param buf Pointer to data buffer.
/// @param nor Number of rows to write. May be zero.
/// @param noc Row size (number of columns). May be zero.

void ffWriteRows(MtxFile_t* file, PTR buf, uint32_t nor, uint32_t noc)
{
   if (noc == 0 || nor == 0) {
      return;
   }
   int ok = 1;
   const size_t rowSizeUsed = ffRowSizeUsed(noc);
   const size_t rowSize = ffRowSize(noc);
   if (rowSizeUsed == rowSize) {
      ok = fwrite(buf, rowSizeUsed, nor, file->file) == nor;
   }
   else {
      const uint8_t* const bufEnd = (const uint8_t*)buf + nor * rowSize;
      for (const uint8_t* b = (const uint8_t*)buf; ok && b < bufEnd; b += rowSize) {
         ok = fwrite(b, rowSizeUsed, 1, file->file) == 1;
      }
   }
   if (!ok) {
      mtxAbort(MTX_HERE, "Cannot write to %s: %s", file->name, strerror(errno));
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads row vectors from a file.
/// @param f Pointer to the file.
/// @param buf Data buffer.
/// @param nor Number of rows to read.
/// @param noc Row size (number of columns).

void ffReadRows(MtxFile_t *f, PTR buf, uint32_t nor, uint32_t noc)
{
   mfValidate(MTX_HERE, f);

   // Handle empty data.
   if (nor == 0 || noc == 0)
      return;

   ffSetField(f->header[0]); // TODO: remove?
   const size_t rowSizeUsed = ffRowSizeUsed(noc);
   const size_t rowSize = ffRowSize(noc);

   // Read rows.
   uint8_t *b = (uint8_t *) buf;
   for (uint32_t i = nor; i > 0; --i) {
      if (fread(b, rowSizeUsed, 1, f->file) != 1) {
         mtxAbort(MTX_HERE,"%s: read error: %s", f->name, strerror(errno));
      }
      b += rowSize;
   }
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
