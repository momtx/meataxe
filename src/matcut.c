////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Cut a rectangular piece out of a matrix.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Extracts a rectangular part of a matrix.
///
/// This function creates a new matrix containing a copy of a rectangular region of the source
/// matrix. The region must not exceed the matrix boundaries.
///
/// See also @ref matDupRows.
///
/// @param src Pointer to the source matrix.
/// @param row0 First row in region.
/// @param col0 First column in region.
/// @param nrows Number of rows to extract.
/// @param ncols Number of columns to extract.
/// @return Pointer to a new matrix containing the specified region, or NULL on error.

Matrix_t* matDupRegion(
   const Matrix_t* src, uint32_t row0, uint32_t col0, uint32_t nrows, uint32_t ncols)
{
   matValidate(MTX_HERE, src);
   if (row0 + nrows > src->nor) {
      mtxAbort(MTX_HERE, "Source row index out of bounds");
   }
   if (col0 + ncols > src->noc) {
      mtxAbort(MTX_HERE, "Source column index out of bounds");
   }

   Matrix_t* result = matAlloc(src->field, nrows, ncols);
   if (nrows == 0 || ncols == 0) {
      return result;
   }

   // Copy the requested region
   PTR s = matGetPtr(src, row0);
   PTR d = result->data;
   for (uint32_t n = nrows; n > 0; --n) {
      if (col0 == 0) {
         ffCopyRow(d, s, ncols);
      }
      else {
         for (uint32_t k = 0; k < ncols; ++k) {
            ffInsert(d, k, ffExtract(s, col0 + k));
         }
      }
      ffStepPtr(&d, ncols);
      ffStepPtr(&s, src->noc);
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Copy a range of rows of a matrix.
///
/// This function creates a new matrix containing a range of consecutive rows of the source matrix.
/// The range must now exceed the matrix's dimensions.
///
/// See also @ref matDupRegion
///
/// @param src Pointer to the source matrix.
/// @param row0 First row in region.
/// @param nrows Number of rows to cut.
/// @return A new matrix containing the specified rows of @em src, or 0 on error.

Matrix_t *matDupRows(const Matrix_t *src, uint32_t row0, uint32_t nrows)
{
   return matDupRegion(src, row0, 0, nrows, src->noc);
}


/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
