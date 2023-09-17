////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Basic matrix functions
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

/// @defgroup mat Matrices over Finite Fields
/// @{
/// @details
/// In the MeatAxe, a matrix over a finite field is represented by a Matrix_t structure.
/// Matrices can be created in many ways, for example
/// - by calling matAlloc(),
/// - by making a copy of an existing matrix with matDup(), or
/// - by reading a matrix from a data file with matRead() or matLoad().
///
/// Matrices that no longer needed must be deleted by calling matFree(). Matrices can consume
/// large amounts of memory, so it always a good idea to delete a matrix as early as possible.
///
/// As with row vectors, row and column indexs are zero-based. For example, in a 3 by 5 matrix
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
      && mat->Magic == MTX_TYPE_MATRIX 
      && mat->Field >= 2 
      && mat->Nor >= 0 
      && mat->Noc >= 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Checks whether the given matrix is valid.
/// If it is not, the program is aborted with an error message.

void matValidate(const struct MtxSourceLocation* src, const Matrix_t *mat)
{
   if (mat == NULL)
      mtxAbort(src,"NULL matrix");
   if (mat->Magic != MTX_TYPE_MATRIX || mat->Field < 2 || mat->Nor < 0 || mat->Noc < 0) {
      mtxAbort(src ? src : MTX_HERE,"Invalid matrix (field=%d, nor=%d, noc=%d)",
            mat->Field, mat->Nor,mat->Noc);
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Create a new matrix.
/// This function creates a new matrix with given dimensions over a given field.
/// @attention
/// To destroy a matrix, use matFree(), not sysFree().
///
/// @param field Field order.
/// @param nor Number of rows.
/// @param noc Number of columns.
/// @return Pointer to the new matrix or NULL on error.

Matrix_t *matAlloc(int field, int nor, int noc)
{
   Matrix_t *m;

   MTX_ASSERT(field >= 2);
   MTX_ASSERT(nor >= 0);
   MTX_ASSERT(noc >= 0);

   // Initialize the data structure
   if (ffSetField(field) != 0) {
      mtxAbort(MTX_HERE,"Cannot select field GF(%d)",field);
      return NULL;
   }

   // Allocate a new Matrix_t structure
   m = ALLOC(Matrix_t);
   if (m == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate Matrix_t structure");
      return NULL;
   }

   // Initialize the data structure
   m->Magic = MTX_TYPE_MATRIX;
   m->Field = field;
   m->Nor = nor;
   m->Noc = noc;
   m->PivotTable = NULL;
   m->Data = ffAlloc(nor, noc);
   if (m->Data == NULL) {
      sysFree(m);
      mtxAbort(MTX_HERE,"Cannot allocate matrix data");
      return NULL;
   }
   return m;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Pointer to a row of a matrix.
/// This function returns a pointer to the specified row of a matrix.
/// Row numbers start from 0. The current row size is not changed.
/// @param mat Pointer to the matrix.
/// @param row Row index.
/// @return Pointer to the selected row or 0 on error.

PTR matGetPtr(const Matrix_t *mat, int row)
{
#ifdef MTX_DEBUG
   matValidate(MTX_HERE, mat);
#ifdef PARANOID
   if ((row < 0) || (row >= mat->Nor))
#else
   if ((row < 0) || (row > mat->Nor + 5))
#endif
   {
      mtxAbort(MTX_HERE,"row=%d: %s",row,MTX_ERR_BADARG);
      return NULL;
   }
#endif
   ffSetField(mat->Field);
   return ffGetPtr(mat->Data, row, mat->Noc);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Delete the pivot table of a matrix.
/// This function deletes the pivot table associated with a matrix. It is used internally,
/// applications should never call this function directly.
/// @param mat Pointer to the matrix.
////

void mat_DeletePivotTable(Matrix_t *mat)
{
   if (mat->PivotTable != NULL) {
      sysFree(mat->PivotTable);
      mat->PivotTable = NULL;
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Delete a matrix.
/// This function frees a matrix which has beed created by matAlloc(). Freeing includes the
/// internal data buffers as well as the Matrix_t structure itself.
/// @param mat Pointer to the matrix.
/// @return 0 on success, -1 on error.

int matFree(Matrix_t *mat)
{
   matValidate(MTX_HERE, mat);
   mat_DeletePivotTable(mat);
   if (mat->Data != NULL) {
      sysFree(mat->Data);
   }
   memset(mat,0,sizeof(Matrix_t));
   sysFree(mat);
   return 0;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
