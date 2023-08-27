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
/// Write row vectors to a file.
/// This function writes |nrows| rows from a data file into a buffer.
/// Unlike |ffWriteRows()|, this function changes the current row size to
/// the appropriate value, which is stored in the |MtxFile_t| object.
/// @param f Pointer to the file.
/// @param buf Data buffer.
/// @param nrows Number of rows to write.
/// @return Number of rows that were actually written. Any value other than count
///    indicates an error.

int mfWriteRows(MtxFile_t *f, PTR buf, int nrows)
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

   // Write rows one by one
   for (i = 0; i < nrows; ++i) {
      if (fwrite(b,ffRowSizeUsed(f->Noc),1,f->File) != 1) {
         break;
      }
      b += ffRowSize(f->Noc);
   }
   if (ferror(f->File)) {
      mtxAbort(MTX_HERE,"%s: Write failed: %S",f->Name);
   }
   return i;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
