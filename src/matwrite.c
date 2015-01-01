////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Write a matrix into a file
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

/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write a matrix to a file.
/// @see MatSave
/// @param mat Pointer to the matrix.
/// @param f Pointer to the file.
/// @return 0 on success, -1 on error.

int MatWrite(const Matrix_t *mat, FILE *f)
{
   long hdr[3];

   if (!MatIsValid(mat)) {
      return -1;
   }
   hdr[0] = mat->Field;
   hdr[1] = mat->Nor;
   hdr[2] = mat->Noc;
   if (SysWriteLong(f,hdr,3) != 3) {
      MTX_ERROR("Cannot write header");
      return -1;
   }
   FfSetField(mat->Field);
   FfSetNoc(mat->Noc);
   if (FfWriteRows(f,mat->Data,mat->Nor) != mat->Nor) {
      return -1;
   }
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write a matrix to a file.
/// This function opens a file, writes a matrix to the file, and closes the
/// file. If a file with the specified name already exists, the old contents
/// of the file are destroyd.
/// To write more than one matrix to a file, use MatWrite().
/// @param mat Pointer to the matrix.
/// @param fn File name.
/// @return 0 on success, -1 on error.

int MatSave(const Matrix_t *mat, const char *fn)
{
   FILE *f;
   int i;

   if (!MatIsValid(mat)) {
      return -1;
   }
   if ((f = SysFopen(fn,FM_CREATE)) == NULL) {
      MTX_ERROR1("Cannot open %s: %S",fn);
      return -1;
   }
   i = MatWrite(mat,f);
   fclose(f);
   if (i != 0) {
      MTX_ERROR1("Cannot write matrix to %s",fn);
   }
   return i;
}


/// @}