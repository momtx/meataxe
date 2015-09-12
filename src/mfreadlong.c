////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Read long integers
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
/// Read long integers from a file.
/// This function reads |count| long integers from a data file into a buffer.
/// If necessary, the data is converted from machine independent format into
/// the format needed by the platform. See |SysReadLong()| for details.
/// @param f Pointer to the file.
/// @param buf Data buffer.
/// @param count Number of integers to read.
/// @return Number of integers that were actually read.
///    Any value other than |count| indicates an error.

int MfReadLong(MtxFile_t *f, long *buf, int count)
{
   int rc;
   if (!MfIsValid(f)) {
      return -1;
   }
   rc = SysReadLong(f->File,buf,count);
   if (rc < 0) {
      MTX_ERROR1("%s: read failed",f->Name);
   }
   return rc;
}


/// @}
