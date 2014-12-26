/* ============================= C MeatAxe ==================================
   File:        $Id: matcut.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Cut a rectangular piece out of a matrix.
   --------------------------------------------------------------------------
   (C) Copyright 1997 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"

   
/* --------------------------------------------------------------------------
   Local data
   -------------------------------------------------------------------------- */

MTX_DEFINE_FILE_INFO


/// @addtogroup mat
/// @{


/// Cut a rectangle out of a matrix.
/// This function creates a new matrix containing a copy of a rectangular region of the
/// source matrix. The region, defined by @em row1, @em col1, @em nrows and @em ncols, 
/// must not exceed the matrix. However, both @em nrows and @em ncols may be -1. In this
/// case the region extends up to the last row or last column, respectivly. For example,
/// to extract the first 10 rows from a matrix independently of the number of columns,
/// you could say
/// @code
/// MatCut(mat,0,0,10,-1)
/// @endcode
/// @see MatCopyRegion MatCutRows
/// @param src Pointer to the matrix.
/// @param row1 First row in region.
/// @param col1 First column in region.
/// @param nrows Number of rows to cut. -1 means as many rows as possible.
/// @param ncols Number of columns to cut. -1 means as many columns as possible.
/// @return Pointer to a new matrix containing the specified region, or 0 on error.

Matrix_t *MatCut(const Matrix_t *src, int row1, int col1, int nrows, int ncols)
{
    Matrix_t *result;
    PTR s, d;
    int n;
    
    /* Check arguments
       --------------- */
#ifdef DEBUG
    if (!MatIsValid(src))
	return NULL;
#endif
    if (nrows == -1)
	nrows = src->Nor - row1;
    if (ncols == -1)
	ncols = src->Noc - col1;
    if (row1 < 0 || nrows < 0 || row1 + nrows > src->Nor)
    {
	MTX_ERROR("Source row index out of bounds");
	return NULL;
    }
    if (col1 < 0 || ncols < 0 || col1 + ncols > src->Noc)
	return MTX_ERROR("Source column index out of bounds"), NULL;

    /* Allocate a new matrix for the result
       ------------------------------------ */
    result = MatAlloc(src->Field,nrows,ncols);
    if (result == NULL)
	return NULL;
    if (nrows == 0)
	return result;

    /* Initialize pointers to the source and destination matrix 
       -------------------------------------------------------- */
    s = MatGetPtr(src,row1);
    d = result->Data;

    /* Copy the requested data
       ----------------------- */
    FfSetNoc(ncols);
    for (n = nrows; n > 0; --n)
    {
	if (col1 == 0)
	    FfCopyRow(d,s);
	else
	{
	    register int k;
	    for (k = 0; k < ncols; ++k)
	    {
#ifdef PARANOID
		FEL f;
		FfSetNoc(src->Noc);
		f = FfExtract(s,col1+k);
		FfSetNoc(ncols);
		FfInsert(d,k,f);
#else
		FfInsert(d,k,FfExtract(s,col1+k));
#endif
	    }
	}
	FfStepPtr(&d);
	s = (PTR)((char *)s + src->RowSize);
    }

    return result;
}




/// Copy a range of rows of a matrix.
/// This function creates a new matrix containing a range of consecutive rows of
/// the source matrix. The range must now exceed the matrix's dimensions. However,
/// @em nrows may be given as -1, meaning "up to the last row".
/// @see MatCopyRegion MatCutRows
/// @param src Pointer to the matrix.
/// @param row1 First row in region.
/// @param nrows Number of rows to cut. -1 means as many rows as possible.
/// @return A new matrix containing the specified rows of @em src, or 0 on error.

Matrix_t *MatCutRows(const Matrix_t *src, int row1, int nrows)
{
    return MatCut(src,row1,0,nrows,-1);
}

/// @}
