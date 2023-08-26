////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Bit Strings, file input
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup bs
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read a bit string from a file.
/// @param f File to read from. Must be open for reading.
/// @return The bit string, 0 on error.

BitString_t *bsRead(FILE *f)
{
   BitString_t *b;
   long hdr[3];             /* File header */
   int size;

   if (sysReadLong32(f,hdr,3) != 3) {
      mtxAbort(MTX_HERE,"Cannot read bit string header");
      return NULL;
   }
   if ((hdr[0] != -3) || (hdr[2] != 0)) {
      mtxAbort(MTX_HERE,"Invalid bit string header (%ld,%ld,%ld)",hdr[0],hdr[1],hdr[2]);
      return NULL;
   }
   size = (int) hdr[1];
   if (size < 0) {
      mtxAbort(MTX_HERE,"Invalid bit string size %d in file)",size);
      return NULL;
   }
   b = bsAlloc(size);
   if (b == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate bit string");
      return NULL;
   }
   if (sysReadLongX(f,b->Data,(b->Size + 7) / 8) != (b->Size + 7) / 8) {
      mtxAbort(MTX_HERE,"Cannot read bit string data");
      bsFree(b);
      return NULL;
   }
   return b;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
