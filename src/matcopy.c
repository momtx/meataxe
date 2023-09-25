////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Copy (part of) a matrix
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
// Contributions by Simon King <simon.king@uni-jena.de>
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Copy a rectangular region of a matrix
/// This function copies a rectangular region of @em src tp @em dest. The source region
/// is defined by its upper left corner and dimensions, the destination region is specified
/// by its upper left corner and has the same dimensions.
/// Both @em nrows and @em ncols can be given as -1. In this case the region extends up to
/// the last row or last column, respectively.
/// The two matrices must be over the same field. Both source and destination region must
/// not exceed the matrices' dimensions. In particular, it is not possible to extend the
/// destination matrix by using %matCopyRegion().
/// @param dest Pointer to the destination matrix.
/// @param destrow Destination row.
/// @param destcol Destination column.
/// @param src Pointer to the source matrix.
/// @param row1 First row in region.
/// @param col1 First column in region.
/// @param nrows Number of rows to copy. -1 means as many rows as possible.
/// @param ncols Number of columns to copy. -1 means as many columns as possible.
/// @return 0 on success, -1 on error.

int matCopyRegion(Matrix_t *dest, int destrow, int destcol,
                  const Matrix_t *src, int row1, int col1, int nrows, int ncols)
{
   PTR s, d;
   int i, k;

   /* Check the arguments
      ------------------- */
   matValidate(MTX_HERE, src);
   matValidate(MTX_HERE, dest);
   if (src->field != dest->field) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
      return -1;
   }
   if (nrows == -1) {
      nrows = src->nor - row1;
   }
   if (ncols == -1) {
      ncols = src->noc - col1;
   }
   if ((row1 < 0) || (nrows < 0) || (row1 + nrows > src->nor)) {
      mtxAbort(MTX_HERE,"Source row index out of range");
      return -1;
   }
   if ((col1 < 0) || (ncols < 0) || (col1 + ncols > src->noc)) {
      mtxAbort(MTX_HERE,"Source column index out of range");
      return -1;
   }
   if ((destrow < 0) || (destrow + nrows > dest->nor)) {
      mtxAbort(MTX_HERE,"Destination row index out of range");
      return -1;
   }
   if ((destcol < 0) || (destcol + ncols > dest->noc)) {
      mtxAbort(MTX_HERE,"Destination column index out of range");
      return -1;
   }

   /* Initialize data pointers
      ------------------------ */
   ffSetField(src->field);
#ifdef MTX_DEBUG
   s = row1 < src->nor ? matGetPtr(src,row1) : NULL;
   d = destrow < dest->nor ? matGetPtr(dest,destrow) : NULL;
#else
   s = matGetPtr(src,row1);
   d = matGetPtr(dest,destrow);
#endif

   /* Copy the rectangle
      ------------------ */
   for (i = row1; i < row1 + nrows; ++i) {
      for (k = col1; k < col1 + ncols; ++k) {
#ifdef MTX_DEBUG
         FEL f;
         f = ffExtract(s,k);
         ffInsert(d,destcol + k - col1,f);
#else
         ffInsert(d,destcol + k - col1,ffExtract(s,k));
#endif
      }
      ffStepPtr(&s, src->noc);
      ffStepPtr(&d, dest->noc);
   }

   mat_DeletePivotTable(dest);

   return 0;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
