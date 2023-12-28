////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Basic matrix functions
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

/// @defgroup mat Matrices over Finite Fields
/// @{
/// In the MeatAxe, a matrix over a finite field is represented by a Matrix_t structure.
/// Matrices can be created in many ways, for example
/// - by calling matAlloc(),
/// - by making a copy of an existing matrix with matDup(), or
/// - by reading a matrix from a data file with matRead() or matLoad().
///
/// Matrices that no longer needed must be deleted by calling matFree(). Matrices can consume
/// large amounts of memory, so it always a good idea to delete a matrix as early as possible.
///
/// Row and column indexes are zero-based. For example, in a 3 by 5 matrix
/// the row index runs from 0 to 2 and the column index runs from 0 to 4.
///
/// A matrix A with entries (a<sub>ij</sub>) is said to be in <b>echelon form</b>
/// if the following conditions are satisfied:
///
/// - Each row has a first non-zero element, called the <b>pivot element</b>.
///   The pivot element may have any value except zero.
/// - If a<sub>ij</sub> is the pivot element of the i-th row, all elements below are zero,
///   i.e., a<sub>ik</sub>=0 for all k>i.
///
/// If a matrix is in echelon form, the column indexes of its pivot elements are
/// called the <b>pivot columns</b> of the matrix. The list of all pivot columns
/// is called the <b>pivot table</b> of the matrix.
/// For a matrix in echelon form the number of rows, @c Nor, is always less than or equal to
/// the number of columns, @c Noc.
/// When different from 0, the @c PivotTable field is a pointer to an array of integers
/// containing first the pivot columns
/// and then the non-pivot columns. This means, the size of the array is
/// always @c Noc: the first @c Nor elements contain the pivot columns
/// and the remaining @c Noc-@c Nor elements contain the non-pivot columns in
/// arbitrary order.

/// @class Matrix_t
/// @brief
/// The @c Data field is a pointer to the matrix elements, organized as an array
/// of rows. Note that the size of a row is different from the number
/// of columns because more than one mark can be packed into a byte, and rows are always
/// padded to a multiple of sizeof(long).
/// Both Nor and Noc may be zero. In this case, Data ist still a valid pointer, but the
/// memory block it points to has size zero. The pivot table is optional and can be NULL.
/// Besides the marks, each matrix carries the field order, the number of rows and the number
/// of columns. There is no global field order or row length as at the kernel layer.

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns true if the given matrix is valid.

int matIsValid(const Matrix_t *mat)
{
   return mat != NULL
      && mat->typeId == MTX_TYPE_MATRIX 
      && mat->field >= 2 
      && mat->nor >= 0 
      && mat->noc >= 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Aborts the program if the passed matrix is not valid.

void matValidate(const struct MtxSourceLocation* src, const Matrix_t* mat)
{
   if (mat == NULL) {
      mtxAbort(src, "NULL matrix");
   }
   if (mat->typeId != MTX_TYPE_MATRIX || mat->field < 2 || mat->nor < 0 || mat->noc < 0) {
      mtxAbort(src ? src : MTX_HERE, "Invalid matrix (field=%d, nor=%lu, noc=%lu)",
         mat->field, (unsigned long)mat->nor, (unsigned long) mat->noc);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Creates a new matrix.
///
/// @param field Field order.
/// @param nor Number of rows.
/// @param noc Number of columns.
/// @return Pointer to the new matrix

Matrix_t *matAlloc(int field, uint32_t nor, uint32_t noc)
{
   Matrix_t *m;

   MTX_ASSERT(field >= 2);
   ffSetField(field);
   m = (Matrix_t*) mmAlloc(MTX_TYPE_MATRIX, sizeof(Matrix_t));
   m->field = field;
   m->nor = nor;
   m->noc = noc;
   m->pivotTable = NULL;
   m->data = ffAlloc(nor, noc);
   return m;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Duplicate a matrix
/// This function creates an independent copy of a matrix.

Matrix_t *matDup(const Matrix_t *src)
{
   matValidate(MTX_HERE, src);
   Matrix_t* m = matAlloc(src->field,src->nor,src->noc);
   memcpy(m->data,src->data,ffSize(src->nor, src->noc));
   // copy pivot table, too?
   return m;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Pointer to a row of a matrix.
/// This function returns a pointer to the specified row of a matrix.
/// Row numbers start from 0. The current row size is not changed.
/// @param mat Pointer to the matrix.
/// @param row Row index.
/// @return Pointer to the selected row or 0 on error.

PTR matGetPtr(const Matrix_t *mat, uint32_t row)
{
#ifdef MTX_DEBUG
   matValidate(MTX_HERE, mat);
#ifdef PARANOID
   if ((row >= mat->nor))
#else
   if (row > mat->nor)
#endif
   {
      mtxAbort(MTX_HERE,"row=%d: %s",row,MTX_ERR_BADARG);
      return NULL;
   }
#endif
   ffSetField(mat->field);
   return ffGetPtr(mat->data, row, mat->noc);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Delete the pivot table of a matrix.
/// This function deletes the pivot table associated with a matrix. It is used internally,
/// applications should never call this function directly.

void mat_DeletePivotTable(Matrix_t *mat)
{
   if (mat->pivotTable != NULL) {
      sysFree(mat->pivotTable);
      mat->pivotTable = NULL;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Deletes a matrix and releases all associated resources.

void matFree(Matrix_t *mat)
{
   matValidate(MTX_HERE, mat);
   mat_DeletePivotTable(mat);
   sysFree(mat->data);
   mat->data = NULL;
   mmFree(mat, MTX_TYPE_MATRIX);
}

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

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Copy a rectangular region of a matrix
/// This function copies a rectangular region of @em src tp @em dest. The source region
/// is defined by its upper left corner and dimensions, the destination region is specified
/// by its upper left corner and has the same dimensions.
/// The two matrices must be over the same field. Both source and destination region must
/// not exceed the matrices' dimensions. In particular, it is not possible to extend the
/// destination matrix by using %matCopyRegion().
//
/// @param dest Pointer to the destination matrix.
/// @param drow Destination row.
/// @param dcol Destination column.
/// @param src Pointer to the source matrix.
/// @param srow First row in region.
/// @param scol First column in region.
/// @param snor Number of rows to copy.
/// @param snoc Number of columns to copy.

void matCopyRegion(
      Matrix_t *dest, uint32_t drow, uint32_t dcol,
      const Matrix_t *src, uint32_t srow, uint32_t scol, uint32_t snor, uint32_t snoc)
{
   // Check the arguments
   matValidate(MTX_HERE, src);
   matValidate(MTX_HERE, dest);
   if (src->field != dest->field) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
   }
   if (srow + snor > src->nor) {
      mtxAbort(MTX_HERE,"Source row index out of range");
   }
   if (scol + snoc > src->noc) {
      mtxAbort(MTX_HERE,"Source column index out of range");
   }
   if (drow + snor > dest->nor) {
      mtxAbort(MTX_HERE,"Destination row index out of range");
   }
   if (dcol + snoc > dest->noc) {
      mtxAbort(MTX_HERE,"Destination column index out of range");
   }

   // Initialize data pointers
   ffSetField(src->field);
#ifdef MTX_DEBUG
   PTR s = srow < src->nor ? matGetPtr(src,srow) : NULL;
   PTR d = drow < dest->nor ? matGetPtr(dest,drow) : NULL;
#else
   PTR s = matGetPtr(src,srow);
   PTR d = matGetPtr(dest,drow);
#endif

   // Copy the rectangle
   for (uint32_t i = srow; i < srow + snor; ++i) {
      for (uint32_t k = scol; k < scol + snoc; ++k) {
#ifdef MTX_DEBUG
         FEL f;
         f = ffExtract(s,k);
         ffInsert(d,dcol + k - scol,f);
#else
         ffInsert(d,dcol + k - scol,ffExtract(s,k));
#endif
      }
      ffStepPtr(&s, src->noc);
      ffStepPtr(&d, dest->noc);
   }

   mat_DeletePivotTable(dest);
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
