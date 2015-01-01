////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Bit Strings, file output
//
// (C) Copyright 1998-2014 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

/// @addtogroup bs
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write a bit string to a file.
/// @param bs The bit string.
/// @param f File to write to. Must be open for writing.
/// @return 0 on success, -1 on error.

int BsWrite(BitString_t *bs, FILE *f)
{
   long file_header[3];

   // check arguments
   if (!BsIsValid(bs)) {
      MTX_ERROR1("bs: %E",MTX_ERR_BADARG);
      return -1;
   }
   if (f == NULL) {
      MTX_ERROR1("f: %E",MTX_ERR_BADARG);
      return -1;
   }

   // write header
   file_header[0] = -3;
   file_header[1] = bs->Size;
   file_header[2] = 0;
   if (SysWriteLong32(f,file_header,3) != 3) {
      MTX_ERROR("Cannot write bit string header");
      return -1;
   }

   // write data
   if (SysWriteLongX(f,bs->Data,(bs->Size + 7) / 8) != (bs->Size + 7) / 8) {
      MTX_ERROR("Cannot write bit string data");
      return -1;
   }

   return 0;
}


/// @}
