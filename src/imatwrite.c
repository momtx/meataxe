////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Write an integer matrix into a file
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

/// @addtogroup imat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write an integer matrix to a file.
/// @see ImatSave
/// @param mat Pointer to the matrix.
/// @param f Pointer to the file.
/// @return 0 on success, -1 on error.

int ImatWrite(const IntMatrix_t *mat, FILE *f)
{
   long hdr[3];

   if (!ImatIsValid(mat)) {
      return -1;
   }
   hdr[0] = -8;     // HACK: T_IMAT in 2.3
   hdr[1] = mat->Nor;
   hdr[2] = mat->Noc;
   if (SysWriteLong(f,hdr,3) != 3) {
      MTX_ERROR("Cannot write header");
      return -1;
   }
   if (SysWriteLong(f,mat->Data,mat->Nor * mat->Noc) != mat->Nor * mat->Noc) {
      return -1;
   }
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write an integer matrix to a file.
/// This function writes an integer matrix to a named file. If the file
/// exists, it is destroyed.
/// @see ImatWrite()
/// @param mat Pointer to the matrix.
/// @param file_name File name.
/// @return 0 on success, -1 on error.

int ImatSave(const IntMatrix_t *mat, const char *file_name)
{
   FILE *f;
   int rc;

   if (!ImatIsValid(mat)) {
      return -1;
   }
   f = SysFopen(file_name,FM_CREATE);
   if (f == NULL) {
      return -1;
   }
   rc = ImatWrite(mat,f);
   fclose(f);
   return rc;
}


/// @}
