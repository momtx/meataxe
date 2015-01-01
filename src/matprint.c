////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Clean a row vector
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

/*
   MTX_DEFINE_FILE_INFO
 */

/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Print a matrix on stdout.
/// This function prints a matrix on the standard output in readable form.
/// If @a name is not 0, the name followed by an equal sign is printed before the matrix.
/// @param name Name to print before the matrix, or 0.
/// @param m Pointer to the matrix.

void MatPrint(const char *name, const Matrix_t *m)
{
   PTR x;
   long i, k;

   if (!MatIsValid(m)) {
      return;
   }
   FfSetField(m->Field);
   FfSetNoc(m->Noc);
   x = m->Data;
   if (name != NULL) { printf("%s=\n",name); }
   for (i = 0; i < m->Nor; ++i) {
      for (k = 0; k < m->Noc; ++k) {
         printf("%d",FfExtract(x,k));
      }
      printf("\n");
      FfStepPtr(&x);
   }
}


/// @}