////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Compare matrices
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


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
/// - Otherwise, the relation is determined by the return value of ffCmpRow() on the first row
///   that is not equal in both matrices.
///
/// In case an error occurs, the return value is -2.
///
/// @param a First matrix.
/// @param b Second matrix.
/// @return 0 if the matrices are equal, ±1 otherwise, -2 on error

int matCompare(const Matrix_t *a, const Matrix_t *b)
{
   int i;

   // check arguments
   matValidate(MTX_HERE,a);
   matValidate(MTX_HERE, b);

   // compare fields and dimensions
   if (a->field > b->field) return 1;
   if (a->field < b->field) return -1;
   if (a->noc > b->noc) return 1;
   if (a->noc < b->noc) return -1;
   if (a->nor > b->nor) return 1;
   if (a->nor < b->nor) return -1;

   // Compare the marks row by row. We cannot use memcmp() on the whole matrix
   // because we must ignore padding bytes.
   ffSetField(a->field);
   for (i = 0; i < a->nor; ++i) {
      PTR pa = matGetPtr(a,i);
      PTR pb = matGetPtr(b,i);
      const int diff = ffCmpRows(pa,pb, a->noc);
      if (diff > 0) return 1;
      if (diff < 0) return -1;
   }

   return 0;	// equal
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
