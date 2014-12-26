/* ============================= C MeatAxe ==================================
   File:        $Id: matorder.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Order of a matrix
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include <string.h>

/* --------------------------------------------------------------------------
   Local data
   -------------------------------------------------------------------------- */

MTX_DEFINE_FILE_INFO 



/// @addtogroup mat
/// @{

/// Order of a matrix.
/// This function calculates the order of a matrix. @em mat must be a 
/// non-singular, square matrix. 
/// Even if @em mat is non-singular, the function may fail. This happens if
/// the order is greater than 1000000, or if the order on any cyclic
/// subspace is greater than 1000.
/// @param mat Pointer to the matrix.
/// @return The order of @em mat, or 1 on error.

int MatOrder(const Matrix_t *mat)
{
    PTR m1, v1, v2, v3;
    PTR basis, bend, bptr;
    int *piv;
    char *done;
    int nor, dim;
    int j1, i;
    FEL f1;
    int tord;
    int ord;
    int flag;

    /* Check the matrix
       ---------------- */
    if (!MatIsValid(mat))
	return -1;
    if (mat->Nor != mat->Noc)
	return MTX_ERROR1("%E",MTX_ERR_NOTSQUARE), -1;

    FfSetField(mat->Field);
    FfSetNoc(mat->Noc);
    nor = mat->Nor;
    m1 = FfAlloc(nor);
    memcpy(m1,mat->Data,FfCurrentRowSize * nor);
    bend = basis = FfAlloc(nor+1);

    piv = NALLOC(int,nor+1);
    done = NALLOC(char,nor);
    memset(done,0,(size_t)nor);
    v1 = FfAlloc(1);
    v2 = FfAlloc(1);
    v3 = FfAlloc(1);
    tord = ord = 1;
    dim = 0;
    j1 = 1;
    while (dim < nor && tord <= 1000 && ord <= 1000000)
    {
	/* Get next start vector
	   --------------------- */
	for (j1 = 0; j1 < nor && done[j1]; ++j1);
	if (j1 >= nor) break;	/* Done! */
	FfMulRow(v1,FF_ZERO);
	FfInsert(v1,j1,FF_ONE);

	/* Calculate order on cyclic subspace
	   ---------------------------------- */
	tord = 0;
	flag = 1;
	FfCopyRow(v3,v1);
	do
	{   FfCopyRow(v2,v3);
	    if (flag)
	    {   FfCopyRow(bend,v3);
		bptr = basis;
		for (i = 0; i < dim; ++i)
		{   
		    f1 = FfExtract(bend,piv[i]);
		    if (f1 != 0)
		    {
			FfAddMulRow(bend,bptr,FfNeg(FfDiv(f1,
			    FfExtract(bptr,piv[i]))));
		    }
		    FfStepPtr(&bptr);
		}
		if ((piv[dim] = FfFindPivot(bend,&f1)) >= 0)
		{   
		    done[piv[dim]] = 1;
		    ++dim;
		    FfStepPtr(&bend);
		}
		else
		    flag = 0;
	    }
	    FfMapRow(v2,m1,nor,v3);
	    ++tord;
	}
	while (tord <= 1000 && FfCmpRows(v3,v1) != 0);

	/* Order = lcm(orders on cyclic subspaces)
	   --------------------------------------- */
	ord = lcm(ord,tord);
    }

    /* Clean up
       -------- */
    SysFree(done);
    SysFree(v1);
    SysFree(v2);
    SysFree(v3);
    SysFree(m1);
    SysFree(basis);
    SysFree(piv);

    if (tord > 1000 || ord > 1000000)
	return -1;
    return ord;
}


/// @}
