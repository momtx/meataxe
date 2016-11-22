////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Compare matrices
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Compare two matrices
/// If the matrices are equal, the return value is 0. Otherwise the return value is positive,
/// if @em a is "greater" than @em b and negative, if @em a is "less" than @em b. The ordering
/// matrices is defined as follows:
///
/// - If the matrices are over different fields, the matrix over the smaller field is smaller.
/// - Otherwise, if the matrices have different number of columns, the matrix with the smaller
///   number of columns is smaller.
/// - Otherwise, if the matrices have different number of rows, the matrix with the smaller
///   number of rows is smaller.
/// - Otherwise, the relation is determined by the return value of FfCmpRow() on the first row
///   that is not equal in both matrices.
///
/// In case an error occurs, the return value is -2.
///
/// @param a First matrix.
/// @param b Second matrix.
/// @return 0 if the matrices are equal, ±1 otherwise, -2 on error

int MatCompare(const Matrix_t *a, const Matrix_t *b)
{
   int i;

   // check arguments
   if (!MatIsValid(a) || !MatIsValid(b)) {
      MTX_ERROR1("%E",MTX_ERR_BADARG);
      return -2;
   }

   // compare fields and dimensions
   if (a->Field > b->Field) return 1;
   if (a->Field < b->Field) return -1;
   if (a->Noc > b->Noc) return 1;
   if (a->Noc < b->Noc) return -1;
   if (a->Nor > b->Nor) return 1;
   if (a->Nor < b->Nor) return -1;

   // Compare the marks row by row. We cannot use memcmp() on the whole matrix
   // because we must ignore padding bytes.
   FfSetField(a->Field);
   FfSetNoc(a->Noc);
   for (i = 0; i < a->Nor; ++i) {
      PTR pa = MatGetPtr(a,i);
      PTR pb = MatGetPtr(b,i);
      const int diff = FfCmpRows(pa,pb);
      if (diff > 0) return 1;
      if (diff < 0) return -1;
   }

   return 0;	// equal
}

/// @}
