////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Bit Strings, file output
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup bs
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write a bit string to a file.
/// @param bs The bit string.
/// @param f File to write to. Must be open for writing.
/// @return 0 on success, -1 on error.

int bsWrite(BitString_t *bs, FILE *f)
{
   long file_header[3];

   // check arguments
   bsValidate(MTX_HERE, bs);
   if (f == NULL) {
      mtxAbort(MTX_HERE,"f: %s",MTX_ERR_BADARG);
      return -1;
   }

   // write header
   file_header[0] = -3;
   file_header[1] = bs->Size;
   file_header[2] = 0;
   if (sysWriteLong32(f,file_header,3) != 3) {
      mtxAbort(MTX_HERE,"Cannot write bit string header");
      return -1;
   }

   // write data
   if (sysWriteLongX(f,bs->Data,(bs->Size + 7) / 8) != (bs->Size + 7) / 8) {
      mtxAbort(MTX_HERE,"Cannot write bit string data");
      return -1;
   }

   return 0;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
