/* ============================= C MeatAxe ==================================
   File:        $Id: matnull.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Matrix null space.
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include <stdlib.h>

MTX_DEFINE_FILE_INFO


/**
!description
    This function calculates the null-space of a matrix. The matrix is passed 
    as first argument, |nor| is the number of rows, |piv| must be a pointer 
    to an array of at least $|nor|+1$ integers, and |nsp| must be a pointer to 
    a square matrix of size |nor|.

    If the function is successfull (non-negative return value),
    - |matrix| is reduced to echelon form,
    - |nsp| contains the null-space in echelon form, and
    - |piv| contains a pivot table for the null space.
    If |flags| is nonzero, the null-space is not reduced to echelon form,
    and the contents of |piv| are undefined.
 ** @see 
 **/

static long znullsp(PTR matrix, long nor, int *piv, PTR nsp, int flags)
{
    PTR x, y, a, b;
    int i;
    long noc = FfNoc;
    long dim;
    FEL f;

    /* Make the identity matrix in <nsp>.
       ---------------------------------- */
    FfSetNoc(nor);
    x = nsp;
    for (i = 0; i < nor; ++i)
    {
	piv[i] = -1;
	FfMulRow(x,FF_ZERO);
	FfInsert(x,i,FF_ONE);
	FfStepPtr(&x);
    }

    /* Gaussian elimination
       -------------------- */
    x = matrix;
    y = nsp;
    for (i = 0; i < nor; ++i)
    {
	PTR xx = matrix, yy = nsp;
	long k, p;

	for (k = 0; k < i; ++k)
	{
	    FfSetNoc(noc);
	    if ((p = piv[k]) >= 0 && (f = FfExtract(x,p)) != FF_ZERO)
	    {
		f = FfNeg(FfDiv(f,FfExtract(xx,p)));
		FfSetNoc(noc);
		FfAddMulRow(x,xx,f);
		FfSetNoc(nor);
		FfAddMulRow(y,yy,f);
	    }
	    FfSetNoc(noc);
	    FfStepPtr(&xx);
	    FfSetNoc(nor);
	    FfStepPtr(&yy);
	}
	FfSetNoc(noc);
	piv[i] = p = FfFindPivot(x,&f);
	FfSetNoc(noc);
	FfStepPtr(&x);
	FfSetNoc(nor);
	FfStepPtr(&y);
    }

    /* Step 2: Reduce the null space to echelon form.
       ---------------------------------------------- */
    dim = 0;
    x = y = nsp;
    a = b = matrix;
    for (i = 0; i < nor; ++i)
    {
	if (piv[i] == -1)
	{
	    FfSetNoc(nor);
	    if (y != x) FfCopyRow(y,x);
	    if (flags)
		++dim;
	    else
	    {
		FfCleanRow(y,nsp,dim,piv);
	        piv[dim++] = FfFindPivot(y,&f);
	    }
	    FfStepPtr(&y);
	}
	else
	{
	    FfSetNoc(noc);
	    if (b != a) FfCopyRow(b,a);
	    FfStepPtr(&b);
	}
	FfSetNoc(nor);
	FfStepPtr(&x);
	FfSetNoc(noc);
	FfStepPtr(&a);
    }

    return dim;
}



/**
 ** @addtogroup mat
 ** @{
 **/

/**
 ** Null-space of a matrix
 ** This function calculates the null-space of a matrix. Unlike MatNullSpace(), this function
 ** modifies the orginal matrix, but uses less memory since no temporary workspace is allocated.
 ** The result is in echelon form.
 ** @param mat Pointer to the matrix.
 ** @param flags If nonzero, the null-space is not reduced to echelon form.
 ** @return Pointer to the null-space of |mat|, or |NULL| on error.
 **/

Matrix_t *MatNullSpace_(Matrix_t *mat, int flags)
{
    long dim;
    Matrix_t *nsp;

    /* Check arguments
       --------------- */
    if (!MatIsValid(mat))
	return NULL;

    /* Allocate workspace
       ------------------ */
    nsp = MatAlloc(mat->Field,mat->Nor,mat->Nor);	/* worst case */
    if (nsp == NULL) 
	return NULL;
    nsp->PivotTable = NREALLOC(nsp->PivotTable,int,mat->Nor);

    /* Calculate the null-space
       ------------------------ */
    FfSetNoc(mat->Noc);
    dim = znullsp(mat->Data,mat->Nor,nsp->PivotTable,nsp->Data,flags);
    if (flags)
    {
	SysFree(nsp->PivotTable);
	nsp->PivotTable = NULL;
    }

    /* Resize the result buffer to its actual size
       ------------------------------------------- */
    nsp->Nor = dim;
    nsp->Data = (PTR) SysRealloc(nsp->Data,nsp->RowSize*dim);

    return nsp;
}



/**
 ** Null-space of a matrix
 ** This function calculates the null-space of a matrix. Unlike MatNullSpace_() and
 ** MatNullSpace__(), this function does not change the original matrix, but it allocates
 ** a temporary copy of the matrix and thus needs more memory.
 ** @param mat Pointer to the matrix.
 ** @return Pointer to the null-space, or 0 on error.
 **/

Matrix_t *MatNullSpace(const Matrix_t *mat)
{
    Matrix_t *tmp, *nsp;

    /* Check arguments
       --------------- */
#ifdef DEBUG
    if (!MatIsValid(mat))
	return NULL;
#endif

    /* Non-destructive null-space
       -------------------------- */
    if ((tmp = MatDup(mat)) == NULL)
	return MTX_ERROR("Cannot duplicate matrix"), NULL;
    nsp = MatNullSpace_(tmp,0);
    MatFree(tmp);
    return nsp;
}


/**
 ** Null-space of a matrix
 ** This function calculates the null-space of a matrix and deletes the original matrix.
 ** @see MatNullSpace_(), MatNullSpace()
 ** @param mat Pointer to the matrix.
 ** @return Pointer to the null-space, or 0 on error.
 **/

Matrix_t *MatNullSpace__(Matrix_t *mat)
{
    Matrix_t *nsp;
    nsp = MatNullSpace_(mat,0);
    MatFree(mat);
    return nsp;
}



/**
 ** @}
 **/
