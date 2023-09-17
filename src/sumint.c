////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Sum and intersection of vector spaces (Zassenhaus algorithm)
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"
#include <string.h>


/// @addtogroup ff
/// @{


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Sum and Intersection of Two Vector Spaces.
/// Given two vector spaces V,W∊Fⁿ, this function calculates the sum and the
/// intersection of the spaces, using the Zassenhaus algorithm. Each of the two spaces
/// is given by a set of generating vectors, which need not be linearly independent.
/// Before calling %SumAndIntersection() the caller must allocate and initialize two
/// workspaces and a pivot table:
/// - Both workspaces must have n₁+n₂ rows, where n₁ and n₂ are the number of generating
///   vectors for the two subspaces.
/// - Workspace 1 must contain the concatenation of the generating sets for the two
///   subspaces. Work space 2 need not be initialized.
/// - The pivot table, must be large enough for at least n₁+n₂ entries. It need not be initialized.
///
/// The variables pointed to by @a nor1 and @a nor2 must contain the numbers n₁ and n₂,
/// respectively. On return, *@a nor1 contains the dimension of V+W, and *@a nor2 contains
/// the dimension of V∩W. The first dim(V+W) rows of @a wrk1 contain a basis of V+W,
/// and a basis of V∩W can be found in @a wrk2 starting at position dim(V+W).
/// Both bases are in echelon form, and @a piv contains the pivot table for the bases.
/// @param noc  Number of columns in «wrk1» and «wrk2»
/// @param wrk1 Workspace 1.
/// @param nor1 Input: number of generators for V, output: dim(V+W).
/// @param nor2 Input: number of generators for W, output: dim(V∩W).
/// @param wrk2 Workspace 2.
/// @param piv Pivot table.
/// @return 0 on success, -1 on error.

int ffSumAndIntersection(int noc, PTR wrk1, uint32_t *nor1, uint32_t *nor2, PTR wrk2, uint32_t *piv)
{
    uint32_t dim1 = *nor1, dim2 = *nor2;
    uint32_t i, k, sumdim;
    PTR x1, x2, y1, y2, sec;


    // Check the arguments.
    if (wrk1 == NULL || nor1 == NULL || nor2== NULL
	|| wrk2 == NULL || piv == NULL)
    {
	mtxAbort(MTX_HERE,"%s",MTX_ERR_BADARG);
	return -1;
    }

    // Set up the workspace 2. Initially, it contains a copy of V1
    for (x2 = wrk2, i = 0; i < dim1 + dim2; ++i, ffStepPtr(&x2, noc)) {
      ffMulRow(x2,FF_ZERO, noc);
    }
    memcpy(wrk2,wrk1,ffSize(dim1, noc));

    // Step 1: Echelonize workspace 1, repeating all operations on workspace 2.
    x1 = y1 = wrk1;
    x2 = y2 = wrk2;
    k = 0;
    for (i = 0; i < dim1 + dim2; ++i, ffStepPtr(&x1, noc), ffStepPtr(&x2, noc))
    {
	FEL f;
	uint32_t p;
	if (ffCleanRowAndRepeat(x1,wrk1,k,noc, piv,x2,wrk2)) {
	   return -1;
	}
	if ((p = ffFindPivot(x1,&f,noc)) == MTX_NVAL)
	   continue;	/* Null row - ignore */
	if (k < i)
	{
	    ffSwapRows(y1,x1, noc);
	    ffSwapRows(y2,x2, noc);
	}
	piv[k++] = p;
	ffStepPtr(&y1, noc);
	ffStepPtr(&y2, noc);
    }
    sumdim = k;	// Dimension of V + W

    // Step 2: Echelonize the basis of the intersection.
    sec = x2 = y2;
    for (i = sumdim; i < dim1 + dim2; ++i, ffStepPtr(&x2, noc))
    {
	FEL f;
	uint32_t p;
	ffCleanRow(x2,sec,k - sumdim,noc,piv + sumdim);
	if ((p = ffFindPivot(x2,&f,noc)) == MTX_NVAL)
	    continue;
	if (i > k)
	    ffCopyRow(y2,x2, noc);
	piv[k++] = p;
	ffStepPtr(&y2, noc);
    }

    *nor1 = sumdim;	    // Dimension of U + W
    *nor2 = k - sumdim;	    // Dimension of the intersection

    return 0;
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
