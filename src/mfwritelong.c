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
/// Write 32-bit integers to a file.
/// This function writes an array of 32-bit integers from @a buffer into a data file.
/// Each integer is written in LSB first format to the file. See @ref sysWrite32 for details.
/// @param f Pointer to the file.
/// @param buf Data buffer.
/// @param count Number of 32-bit integers to write.

void mfWrite32(MtxFile_t *f, const void *buf, int count)
{
   mfValidate(f);
   sysWrite32(f->File,buf,count);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read 32-bit integers from a file.
/// This function reads @a count 32-bit integers from a data file into a buffer.
/// Each integer is converted from file format (little-endian) into naive format.
/// See @ref sysread32 for details.

void mfRead32(MtxFile_t* f, void* buf, int count)
{
   mfValidate(f);
   sysRead32(f->File,buf,count);
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
