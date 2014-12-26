////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Basic matrix functions
//
// (C) Copyright 1998-2014 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

/// @defgroup mat Matrices over Finite Fields
/// @{
/// @details
/// In the MeatAxe, a matrix over a finite field is represented by a Matrix_t structure.
/// Matrices can be created in many ways, for example
/// - by calling MatAlloc(),
/// - by making a copy of an existing matrix with MatDup(), or
/// - by reading a matrix from a data file with MatRead() or MatLoad().
///
/// Matrices that no longer needed must be deleted by calling MatFree(). Matrices can consume
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
/// of rows. Note that the size of a row (@c RowSize field) is different from the number
/// of columns because more than one mark can be packed into a byte, and rows are always
/// padded to a multiple of sizeof(long).
/// Both Nor and Noc may be zero. In this case, Data ist still a valid pointer, but the
/// memory block it points to has size zero. The pivot table is optional and can be NULL.
/// Besides the marks, each matrix carries the field order, the number of rows and the number
/// of columns. There is no global field order or row length as at the kernel layer.

   
/* --------------------------------------------------------------------------
   Local data
   -------------------------------------------------------------------------- */

MTX_DEFINE_FILE_INFO;

#define MAT_MAGIC 0x6233af91

   
   
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Check if the matrix is valid.
/// This function checks if the argument is a pointer to a valid matrix. If the matrix is o.k.,
/// the function returns 1.  Otherwise, an error is signalled and, if the error handler does not
/// terminate the program, the function returns 0.
/// @param mat Matrix to check.
/// @return 1 if the matrix is valid, 0 otherwise.

int MatIsValid(const Matrix_t *mat)
{
    if (mat == NULL)
    {
	MTX_ERROR("NULL matrix");
	return 0;
    }
    if (mat->Magic != MAT_MAGIC || mat->Field < 2 || mat->Nor < 0 || 
	mat->Noc < 0)
    {
	MTX_ERROR3("Invalid matrix (field=%d, nor=%d, noc=%d)",mat->Field,
	    mat->Nor,mat->Noc);
	return 0;
    }
    return 1;

}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Create a new matrix.
/// This function creates a new matrix with given dimensions over a given field.
/// @attention
/// To destroy a matrix, use MatFree(), not SysFree().
/// @param field Field order.
/// @param nor Number of rows.
/// @param noc Number of columns.
/// @return Pointer to the new matrix or 0 on error.

Matrix_t *MatAlloc(int field, int nor, int noc)
{
    Matrix_t *m;

    MTX_VERIFY(field >= 2);
    MTX_VERIFY(nor >= 0);
    MTX_VERIFY(noc >= 0);

    /* Allocate a new Matrix_t structure
       --------------------------------- */
    m = ALLOC(Matrix_t);
    if (m == NULL)
    {
	MTX_ERROR("Cannot allocate Matrix_t structure");
	return NULL;
    }

    /* Initialize the data structure
       ----------------------------- */
    if (FfSetField(field) != 0)
    {
	MTX_ERROR1("Cannot select field GF(%d)",field);
	SysFree(m);
	return NULL;
    }
    FfSetNoc(noc);
    m->Magic = MAT_MAGIC;
    m->Field = field;
    m->Nor = nor;
    m->Noc = noc;
    m->PivotTable = NULL;
    m->Data = FfAlloc(nor);
    m->RowSize = FfCurrentRowSize;
    if (m->Data == NULL)
    {
	SysFree(m);
	MTX_ERROR("Cannot allocate matrix data");
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

PTR MatGetPtr(const Matrix_t *mat, int row)
{
#ifdef DEBUG
    if (!MatIsValid(mat))
	return NULL;
#ifdef PARANOID
    if (row < 0 || row >= mat->Nor)
#else
    if (row < 0 || row > mat->Nor + 5)
#endif
    {
	MTX_ERROR2("row=%d: %E",row,MTX_ERR_BADARG);
	return NULL;
    }
#endif
    return (PTR)((char *) mat->Data + mat->RowSize * row);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Delete the pivot table of a matrix.
/// This function deletes the pivot table associated with a matrix. It is used internally,
/// applications should never call this function directly.
/// @param mat Pointer to the matrix.
////

void Mat_DeletePivotTable(Matrix_t *mat)
{
    if (mat->PivotTable != NULL)
    {
	SysFree(mat->PivotTable);
	mat->PivotTable = NULL;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Delete a matrix.
/// This function frees a matrix which has beed created by MatAlloc(). Freeing includes the
/// internal data buffers as well as the Matrix_t structure itself.
/// @param mat Pointer to the matrix.
/// @return 0 on success, -1 on error.

int MatFree(Matrix_t *mat)
{
    if (!MatIsValid(mat))
	return -1;
    Mat_DeletePivotTable(mat);
    if (mat->Data != NULL)
	SysFree(mat->Data);
    memset(mat,0,sizeof(Matrix_t));
    SysFree(mat);
    return 0;
}


/// @}
