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

MTX_DEFINE_FILE_INFO

/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Clean a matrix.
/// This function "cleans" a matrix with a space, i.e., it adds suitable linear combinations
/// of the rows in @em sub to the rows of @em mat such that all pivot columns in @em mat are
/// zero. Both matrices must be over the same field and have the same number of colums.
/// The second matrix, @em sub, must be in echelon form. The cleaned matrix is reduced to
/// echelon form.
/// @return Rank of the cleaned matrix, or -1 on error.

int MatClean(Matrix_t *mat, const Matrix_t *sub)
{
   int i;

   /* Check the arguments
      ------------------- */
   if (!MatIsValid(mat) || !MatIsValid(sub)) {
      return -1;
   }
   if ((mat->Field != sub->Field) || (mat->Noc != sub->Noc)) {
      MTX_ERROR1("%E",MTX_ERR_INCOMPAT);
      return -1;
   }
   if (sub->PivotTable == NULL) {
      MTX_ERROR1("Subspace: %E",MTX_ERR_NOTECH);
      return -1;
   }

   /* Clean
      ----- */
   FfSetNoc(mat->Noc);
   for (i = 0; i < mat->Nor; ++i) {
      PTR m = MatGetPtr(mat,i);
      FfCleanRow(m,sub->Data,sub->Nor,sub->PivotTable);
   }

   /* Reduce to echelon form
      ---------------------- */
   if (MatEchelonize(mat) < 0) {
      return -1;
   }
   return mat->Nor;
}


/// @}