/* ============================= C MeatAxe ==================================
   File:        $Id: matcopy.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
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

/// Copy a rectangular region of a matrix
/// This function copies a rectangular region of @em src tp @em dest. The source region
/// is defined by its upper left corner and dimensions, the destination region is specified
/// by its upper left corner and has the same dimensions.
/// Both @em nrows and @em ncols can be given as -1. In this case the region extends up to
/// the last row or last column, respectively.
/// The two matrices must be over the same field. Both source and destination region must
/// not exceed the matrices' dimensions. In particular, it is not possible to extend the
/// destination matrix by using %MatCopyRegion().
/// @param dest Pointer to the destination matrix.
/// @param destrow Destination row.
/// @param destcol Destination column.
/// @param src Pointer to the source matrix.
/// @param row1 First row in region.
/// @param col1 First column in region.
/// @param nrows Number of rows to copy. -1 means as many rows as possible.
/// @param ncols Number of columns to copy. -1 means as many columns as possible.
/// @return 0 on success, -1 on error.

int MatCopyRegion(Matrix_t *dest, int destrow, int destcol, 
    const Matrix_t *src, int row1, int col1, int nrows, int ncols)

{
    PTR s, d;
    int i, k;

    /* Check the arguments
       ------------------- */
    if (!MatIsValid(src) || !MatIsValid(dest))
	return -1;
    if (src->Field != dest->Field)
	return MTX_ERROR1("%E",MTX_ERR_INCOMPAT), -1;
    if (nrows == -1)
	nrows = src->Nor - row1;
    if (ncols == -1)
	ncols = src->Noc - col1;
    if (row1 < 0 || nrows < 0 || row1 + nrows > src->Nor)
    {
	MTX_ERROR("Source row index out of range");
	return -1;
    }
    if (col1 < 0 || ncols < 0 || col1 + ncols > src->Noc)
    {
	MTX_ERROR("Source column index out of range");
	return -1;
    }
    if (destrow < 0 || destrow + nrows > dest->Nor)
    {
	MTX_ERROR("Destination row index out of range");
	return -1;
    }
    if (destcol < 0 || destcol + ncols > dest->Noc)
    {
	MTX_ERROR("Destination column index out of range");
	return -1;
    }

    /* Initialize data pointers
       ------------------------ */
    FfSetField(src->Field);
#ifdef PARANOID
    s = row1 < src->Nor ? MatGetPtr(src,row1) : NULL;
	d = destrow < dest->Nor ? MatGetPtr(dest,destrow) : NULL;
#else
    s = MatGetPtr(src,row1);
    d = MatGetPtr(dest,destrow);
#endif

    /* Copy the rectangle
       ------------------ */
    for (i = row1; i < row1 + nrows; ++i)
    {
	for (k = col1; k < col1 + ncols; ++k)
	{
#ifdef PARANOID
	    FEL f;
	    FfSetNoc(src->Noc);
	    f = FfExtract(s,k);
	    FfSetNoc(dest->Noc);
	    FfInsert(d,destcol+k-col1,f);
#else
	    FfInsert(d,destcol+k-col1,FfExtract(s,k));
#endif
	}
	/*s = FfGetPtr(s,1,src->Noc);*/
	s = (PTR)((char *)s + src->RowSize);
	/*d = FfGetPtr(d,1,dest->Noc);*/
	d = (PTR)((char *)d + dest->RowSize);
    }

    Mat_DeletePivotTable(dest);

    return 0;
}

/// @}
