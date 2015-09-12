////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Reduce a matrix to semi echelon form.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include <meataxe.h>
#include <stdlib.h>

MTX_DEFINE_FILE_INFO


/// @addtogroup ff
/// @{


/// Clean Row.
/// This function performs a Gaussian elimination, i.e., it adds suitable multiples of the
/// rows of @em matrix to @em row such that all pivot positions are zero. @em piv is the pivot
/// table for @em matrix. As usual, all indexes are 0-based, i.e., <tt>piv[0]</tt> is the
/// pivot column of the first row, and for a unit matrix we have <tt>piv[0]==0</tt>.
/// The field and row size must have been set before calling this function.
/// @param row The row to be cleaned.
/// @param matrix Pointer to the matrix.
/// @param nor Number of rows of the matrix.
/// @param piv The pivot table.

void FfCleanRow(PTR row, PTR matrix, int nor, const int *piv)
{
    int i;
    PTR x;

    for (i=0, x=matrix; i < nor; ++i, FfStepPtr(&x))
    {
        FEL f = FfExtract(row,piv[i]);
        if (f != FF_ZERO)
	    FfAddMulRow(row,x,FfNeg(FfDiv(f,FfExtract(x,piv[i]))));
    }
}


/// Clean Row and Record Operations.
/// This function works like FfCleanRow(), but it stores a record of the operations performed
/// in @em row2. @em row2 must be a row of at least @em nor entries.  On return, @em row2
/// contains the coefficients by which the rows of @em mat were multiplied and then subtracted
/// from @em row.
/// Before calling %FfCleanRow2(), the caller must initialize @em row2 to zero. Otherwise the
/// results are undefined.
/// @param row Pointer to row to be cleaned.
/// @param mat Matrix to clean with.
/// @param nor Number of rows.
/// @param piv Pivot table for @em matrix.
/// @param row2 Pointer to row where the operations are recorded.
/// @return Always 0.

void FfCleanRow2(PTR row, PTR mat, int nor, const int *piv, PTR row2)
{
    int i;
    PTR x;

    if (row2 == NULL || piv == NULL)
    {
	MTX_ERROR1("%E",MTX_ERR_BADARG);
	return;
    }
    for (i = 0, x = mat; i < nor; ++i, FfStepPtr(&x))
    {
	FEL f = FfExtract(row,piv[i]);
	if (f != FF_ZERO) 
	{
	    f = FfDiv(f,FfExtract(x,piv[i]));
	    FfAddMulRow(row,x,FfNeg(f));
	    FfInsert(row2,i,f);
	}
    }
}



/// Clean Row and Repeat Operations.
/// This function works like FfCleanRow(), but repeats all operations on
/// a second row/matrix. 
/// @param row Pointer to row to be cleaned.
/// @param mat Matrix to clean with.
/// @param nor Number of rows.
/// @param piv Pivot table for @em mat.
/// @param row2 Pointer to the second row to be cleaned.
/// @param mat2 Matrix to the second matrix.
/// @return Always 0.

void FfCleanRowAndRepeat(PTR row, PTR mat, int nor, const int *piv, PTR row2, PTR mat2)
{
    int i;
    PTR x, x2;

#ifdef DEBUG
    if (row2 == NULL || piv == NULL || row2 == NULL || mat2 == NULL)
    {
	MTX_ERROR1("%E",MTX_ERR_BADARG);
	return;
    }
#endif
    for (i = 0, x = mat, x2 = mat2; i < nor; ++i, FfStepPtr(&x), FfStepPtr(&x2))
    {
	FEL f = FfExtract(row,piv[i]);
	if (f != FF_ZERO) 
	{
	    f = FfNeg(FfDiv(f,FfExtract(x,piv[i])));
	    FfAddMulRow(row,x,f);
	    FfAddMulRow(row2,x2,f);
	}
    }
}


/// @}
