////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Bit Strings, file output
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup bs
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Writes a bit string to a file.
/// @param bs The bit string.
/// @param f File to write to. Must be open for writing.

void bsWrite(BitString_t *bs, FILE *file)
{
   // check arguments
   bsValidate(MTX_HERE, bs);
   MTX_ASSERT(file != NULL);

   // write header
   uint32_t fileHeader[3] = {MTX_TYPE_BITSTRING, bs->Size, 0};
   sysWrite32(file, fileHeader, 3);

   // write data
   sysWrite8(file, bs->Data, sysPad(bs->Size, 8));
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads a bit string from a file.
/// @param f File to read from. Must be open for reading.
/// @return The bit string.

BitString_t *bsRead(FILE *file)
{
   MTX_ASSERT(file != NULL);

   uint32_t fileHeader[3];
   sysRead32(file, fileHeader, 3);
   if (fileHeader[0] != MTX_TYPE_BITSTRING || fileHeader[2] != 0) {
      mtxAbort(MTX_HERE,"Invalid bit string header (%lu,%lu,%lu)",
            (unsigned long) fileHeader[0],
            (unsigned long) fileHeader[1],
            (unsigned long) fileHeader[2]);
   }

   BitString_t *bs = bsAlloc(fileHeader[1]);
   sysRead8(file, bs->Data, sysPad(bs->Size, 8));
   return bs;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
