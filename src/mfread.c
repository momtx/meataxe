////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Read rows from a file
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup mf
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read row vectors from a file.
/// This function reads @a nrows rows from a data file into a buffer.
/// Unlike ffReadRows(), this function changes the current row size to
/// the appropriate value, which is stored in the MtxFile_t object.
/// @param f Pointer to the file.
/// @param buf Data buffer.
/// @param nrows Number of rows to read.
/// @return Number of rows that were actually read. Any value other than count indicates an error.

int mfReadRows(MtxFile_t *f, PTR buf, int nrows)
{
   int i;
   register char *b = (char *) buf;

   if (!mfIsValid(f)) {
      return -1;
   }

   // Handle special case Noc=0
   if (f->Noc == 0) {
      return nrows;
   }

   // Read rows one by one
   for (i = 0; i < nrows; ++i) {
      if (fread(b,ffRowSizeUsed(f->Noc),1,f->File) != 1) { break; }
      b += ffRowSize(f->Noc);
   }
   if (ferror(f->File)) {
      mtxAbort(MTX_HERE,"%s: Read failed: %S",f->Name);
   }
   return i;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
