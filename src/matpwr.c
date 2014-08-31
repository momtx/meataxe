/* ============================= C MeatAxe ==================================
   File:        $Id: matpwr.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Power of a matrix.
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include <string.h>

MTX_DEFINE_FILE_INFO

/**
 ** @addtogroup mat
 ** @{
 **/



/* ------------------------------------------------------------------
   matpower() - Calculate the <n>-th power of <mat> using the
	binary method.
   ------------------------------------------------------------------ */

static void matpwr_(long n, PTR inp, PTR out, PTR tmp2)

{
    PTR x, y;
    long i;
    int first = 1;

    while (n > 0)
    {	if (n % 2 == 1)
	{
	    if (first)
	    {
		memcpy(out,inp,FfCurrentRowSize * FfNoc);
		first = 0;
	    }
	    else
	    {
		x = out;
		for (i = 1; i <= FfNoc; ++i)
		{
		    FfMapRow(x,inp,FfNoc,tmp2);
		    FfCopyRow(x,tmp2);
		    FfStepPtr(&x);
		}
	    }
	}
	if (n == 1)
	    break;
	x = inp;
	y = tmp2;
	for (i = 1; i <= FfNoc; ++i)
	{
	    FfMapRow(x,inp,FfNoc,y);
	    FfStepPtr(&x);
	    FfStepPtr(&y);
	}
	memcpy(inp,tmp2,FfCurrentRowSize * FfNoc);
	n /= 2;
    }
}



/**
!section obj.mat
 ** Power of a matrix.
!synopsis 
    Matrix_t *MatPower(const Matrix_t *mat, long n);
 ** @param mat
    Pointer to the matrix.
 ** @param n
    Exponent.
 ** @return
    |n|-th power of |mat|, or |NULL| on error.
!description
    This function calculates the $n$-th power of a matrix, using the binary
    method. This method is generally faster than multiplying the matrix $n$
    times by itself. On the other hand, a third matrix is temporarily created 
    in addition to the original matrix and the result matrix.
    The cases $n=0$ and $n=1$ are handled separately, avoiding unnecessary
    memory allocation and calculation. 

    Negative exponents are not allowed. To calculate a negative power, you 
    must first invert the matrix with |MatInverse()| and then call |MatPower()|
    with the inverted matrix and a positive exponent.
 ** @see MatMul MatInverse
 **/

Matrix_t *MatPower(const Matrix_t *mat, long n)

{
    Matrix_t *result;
    PTR tmp, tmp2;

    /* Check the arguments
       ------------------- */
    if (!MatIsValid(mat))
	return NULL;
    if (mat->Nor != mat->Noc)
    {
	MTX_ERROR1("MatPower(): %E",MTX_ERR_NOTSQUARE);
	return NULL;
    }

    /* Handle special cases n = 0 and n = 1
       ------------------------------------ */
    if (n == 0) 
	return MatId(mat->Field,mat->Nor);
    else if (n == 1)
	return MatDup(mat);

    FfSetField(mat->Field);
    FfSetNoc(mat->Noc);
    tmp = FfAlloc(FfNoc);
    memcpy(tmp,mat->Data,FfCurrentRowSize * FfNoc);
    tmp2 = FfAlloc(FfNoc);
    result = MatAlloc(mat->Field,mat->Nor,mat->Noc);
    if (result != NULL)
	matpwr_(n,tmp,result->Data,tmp2);
    SysFree(tmp);
    SysFree(tmp2);
    return result;
}



/**
 ** @}
 **/
