////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Write long integers
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

/// @addtogroup mf
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write long integers to a file.
/// This function writes |count| long integers from buffer into a data file.
/// If necessary, the data is converted from machine independent format into
/// the format needed by the platform. See |SysWriteLong()| for details.
/// @param f Pointer to the file.
/// @param buf Data buffer.
/// @param count Number of integers to write.
/// @return Number of integers that were actually written. Any value other than count
///    indicates an error.

int MfWriteLong(MtxFile_t *f, const long *buf, int count)
{
   int rc;
   if (!MfIsValid(f)) {
      return -1;
   }
   rc = SysWriteLong(f->File,buf,count);
   if (rc != count) {
      MTX_ERROR1("%s: write failed",f->Name);
   }
   return rc;
}


/// @}
