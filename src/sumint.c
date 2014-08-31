/* ============================= C MeatAxe ==================================
   File:        $Id: sumint.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Sum and intersection of vector spaces (Zassenhaus algorithm)
   --------------------------------------------------------------------------
   (C) Copyright 1999 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include <string.h>

MTX_DEFINE_FILE_INFO

/**
 ** @addtogroup ff2
 ** @{
 **/


/**
 ** Sum and Intersection of Two Vector Spaces.
 ** Given two vector spaces V,W∊F<sup>n</sup>, this function calculates the sum and the
 ** intersection of the spaces, using the Zassenhaus algorithm. Each of the two spaces
 ** is given by a set of generating vectors, which need not be linearly independent.
 ** Before calling %SumAndIntersection() the caller must allocate and initialize two
 ** workspaces and a pivot table:
 ** - Both workspaces must have n₁+n₂ rows, where n₁ and n₂ are the number of generating
 **   vectos for the two subspaces.
 ** - Workspace 1 must contain the concatenation of the generating sets for the two
 **   subspaces. Work space 2 need not be initialized.
 ** - The pivot table, must be large enough for at least n₁+n₂ entries. It need not be initialized.
 **
 ** The variables pointed to by @a nor1 and @a nor2 must contain the numbers n₁ and n₂,
 ** respectively. On return, *@a nor1 contains the dimension of V+W, and *@a nor2 contains
 ** the dimension of V∩W. The first dim(V+W) rows of @a wrk1 contain a basis of V+W,
 ** and a basis of V∩W can be found in @a wrk2 starting at position dim(V+W).
 ** Both bases are in echelon form, and @a piv contains the pivot table for the bases.
 ** @param wrk1 Workspace 1.
 ** @param nor1 Input: number of generators for V, output: dim(V+W).
 ** @param nor2 Input: number of generators for W, output: dim(V∩W).
 ** @param wrk2 Workspace 2.
 ** @param piv Pivot table.
 ** @return 0 on success, -1 on error.
 **/

int FfSumAndIntersection(PTR wrk1, int *nor1, int *nor2, PTR wrk2, int *piv)
{
    int dim1 = *nor1, dim2 = *nor2;
    int i, k, sumdim;
    PTR x1, x2, y1, y2, sec;


    /* Check the arguments.
       -------------------- */
    if (wrk1 == NULL || nor1 == NULL || nor2== NULL
	|| wrk2 == NULL || piv == NULL)
    {
	MTX_ERROR1("%E",MTX_ERR_BADARG);
	return -1;
    }

    /* Set up the workspace 2. Initially, it contains a copy of V1.
       ------------------------------------------------------------ */
    for (x2 = wrk2, i = 0; i < dim1 + dim2; ++i, FfStepPtr(&x2))
	FfMulRow(x2,FF_ZERO);
    memcpy(wrk2,wrk1,dim1 * FfCurrentRowSize);

    /* Step 1: Echelonize workspace 1, repeating 
       all operations on workspace 2.
       ----------------------------------------- */
    x1 = y1 = wrk1;
    x2 = y2 = wrk2;
    k = 0;
    for (i = 0; i < dim1 + dim2; ++i, FfStepPtr(&x1), FfStepPtr(&x2))
    {
	FEL f;
	int p;
	FfCleanRowAndRepeat(x1,wrk1,k,piv,x2,wrk2);
	if ((p = FfFindPivot(x1,&f)) < 0)
	    continue;	/* Null row - ignore */
	if (k < i)
	{
	    FfSwapRows(y1,x1);
	    FfSwapRows(y2,x2);
	}
	piv[k++] = p;
	FfStepPtr(&y1);
	FfStepPtr(&y2);
    }
    sumdim = k;	/* Dimension of V + W */

    /* Step 2: Echelonize the basis of the intersection.
       ------------------------------------------------- */
    sec = x2 = y2;
    for (i = sumdim; i < dim1 + dim2; ++i, FfStepPtr(&x2))
    {
	FEL f;
	int p;
	FfCleanRow(x2,sec,k - sumdim,piv + sumdim);
	if ((p = FfFindPivot(x2,&f)) < 0)
	    continue;
	if (i > k)
	    FfCopyRow(y2,x2);
	piv[k++] = p;
	FfStepPtr(&y2);
    }

    *nor1 = sumdim;	    /* Dimension of U + W */
    *nor2 = k - sumdim;	    /* Dimension of the intersection */

    return 0;
}

/**
 ** @}
 **/
