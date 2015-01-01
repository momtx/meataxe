////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Trace of a matrix
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Trace of a Matrix.
/// This function calculates the sum of all diagonal elements of a matrix.
/// Note that the matrix need not be square.
/// @param mat Pointer to the matrix.
/// @return Trace of @a mat, @c FF_ZERO on error.

FEL MatTrace(const Matrix_t *mat)
{
   int i;
   PTR x;
   int maxi;
   FEL trace = FF_ZERO;

   /* Check the argument
      ------------------ */
#ifdef DEBUG
   if (!MatIsValid(mat)) {
      return FF_ZERO;
   }
#endif

   maxi = mat->Nor > mat->Noc ? mat->Noc : mat->Nor;
   FfSetField(mat->Field);
   FfSetNoc(mat->Noc);
   for (i = 0, x = mat->Data; i < maxi; ++i, FfStepPtr(&x)) {
      trace = FfAdd(trace,FfExtract(x,i));
   }
   return trace;
}


/// @}