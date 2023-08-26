////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Reduce a matrix to semi echelon form.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
// Contributions by Simon King <simon.king@uni-jena.de>
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>


/// @addtogroup ff
/// @{

/// Clean Row.
/// This function performs a Gaussian elimination, i.e., it adds suitable multiples of the
/// rows of @em matrix to @em row such that all pivot positions are zero. @em piv is the pivot
/// table for @em matrix. As usual, all indexes are 0-based, i.e., <tt>piv[0]</tt> is the
/// pivot column of the first row, and for a unit matrix we have <tt>piv[0]==0</tt>.
/// The field must have been set before calling this function.
/// @param row The row to be cleaned.
/// @param matrix Pointer to the matrix.
/// @param nor Number of rows of the matrix.
/// @param nor Number of columns in @a row and @a matrix.
/// @param piv The pivot table.

void ffCleanRow(PTR row, PTR matrix, int nor, int noc, const int *piv)
{
   PTR x;
   int i;

   MTX_ASSERT(ffNoc == noc,);
   for (i = 0, x = matrix; i < nor; ++i, ffStepPtr(&x, noc)) {
      const int pivi = piv[i];
      const FEL f = ffExtract(row,pivi);
      if (f != FF_ZERO) {
         ffAddMulRowPartial(row, x, ffNeg(ffDiv(f, ffExtract(x, pivi))), pivi);
      }
   }
}


/// Clean Row and Record Operations.
/// This function works like ffCleanRow(), but it stores a record of the operations performed
/// in @em row2. @em row2 must be a row of at least @em nor entries.  On return, @em row2
/// contains the coefficients by which the rows of @em mat were multiplied and then subtracted
/// from @em row.
/// Before calling %ffCleanRow2(), the caller must initialize @em row2 to zero. Otherwise the
/// results are undefined.
///
/// @param row Pointer to row to be cleaned.
/// @param mat Matrix to clean with.
/// @param nor Number of rows.
/// @param piv Pivot table for @em matrix.
/// @param row2 Pointer to row where the operations are recorded.
/// @return 0 on success, -1 on error (if any of the pointer arguments is NULL).

int ffCleanRow2(PTR row, PTR mat, int nor, int noc, const int *piv, PTR row2)
{
   int i;
   PTR x;

   if (row == NULL || mat == NULL || row2 == NULL || piv == NULL) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_BADARG);
      return -1;
   }
   for (i = 0, x = mat; i < nor; ++i, ffStepPtr(&x, noc)) {
      FEL f = ffExtract(row,piv[i]);
      if (f != FF_ZERO) {
         f = ffDiv(f,ffExtract(x,piv[i]));
         ffAddMulRow(row,x,ffNeg(f));
         ffInsert(row2,i,f);
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int ffCleanRowAndRepeat(PTR row, PTR mat, int nor, int noc, const int *piv, PTR row2, PTR mat2)
{
   int i;
   PTR x, x2;

   MTX_ASSERT(row != NULL,-1);
   MTX_ASSERT(row2 != NULL,-1);
   MTX_ASSERT(mat2 != NULL,-1);
   MTX_ASSERT(piv != NULL,-1);
   MTX_ASSERT(ffNoc == noc, -1);
   for (i = 0, x = mat, x2 = mat2; i < nor; ++i, ffStepPtr(&x,noc), ffStepPtr(&x2,noc)) {
      FEL f = ffExtract(row,piv[i]);
      if (f != FF_ZERO) {
         f = ffNeg(ffDiv(f,ffExtract(x,piv[i])));
         ffAddMulRow(row,x,f);
         ffAddMulRow(row2,x2,f);
      }
   }
   return 0;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
