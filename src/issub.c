////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Compare vector spaces
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
/// @addtogroup algo
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Vector space incidence relation.
/// This function compares two vector spaces U,V≤F<sup>n</sup> and returns 1
/// if U≤V, or 0 otherwise.
/// The spaces to compare are given as matrices where the rows of @a sub generate U
/// and @a space is a basis of V.
/// Thus, @a sub and @a space must be matrices over the same field and with the same
/// number of columns. @a space must be in full echelon form, but there is no further
/// restriction on @a sub. In particular, the rows of @a sub may be linearly dependent.
///
/// In normal mode, @a ngen is 0. Then, the rows of @a sub are checked one-by-one
/// to see if they are in the vector space generated by @a space. If this test passes
/// for each row of @a sub, the return value is 1, otherwise it is 0. 
///
/// If @a ngen is different from zero, the function assumes that U is  generated by the
/// first @a ngen rows, and only this number of rows are checked. For example, U may be
/// a cyclic submodule, which is generated from a single vector under the action of an
/// algebra. Then, one only needs to check if the generating vector, which 
/// is supposed to be the first row of @a sub, is contained in V.
///
/// @param sub Generating set for the first space, U.
/// @param space Pointer to a basis for the second space, V.
/// @param ngen Number of vectors in @a sub to consider.
/// @return 1 if U≤V, 0 otherwise

int IsSubspace(const Matrix_t *sub, const Matrix_t *space, int ngen)
{
    long nvec, spcdim, i;
    PTR tmp, y;

    matValidate(MTX_HERE, sub);
    matValidate(MTX_HERE, space);
    if (sub->field != space->field || sub->noc != space->noc)
	mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
    if (space->pivotTable == NULL)
	mtxAbort(MTX_HERE,"space: %s",MTX_ERR_NOTECH);

    /* Decide how many vectors from <sub> we shall check
       ------------------------------------------------- */
    ffSetField(space->field);
    spcdim = space->nor;
    nvec = sub->nor;
    if (ngen > 0 && ngen < nvec)
	nvec = ngen;

    tmp = ffAlloc(1, space->noc);

    /* Check if vectors from <sub> are in the span of <space>
       ------------------------------------------------------ */
    for (i = nvec, y = sub->data; i > 0; --i, ffStepPtr(&y, space->noc))
    {
	FEL f;
	ffCopyRow(tmp,y, space->noc);
	ffCleanRow(tmp,space->data,spcdim,space->noc, space->pivotTable);
        if (ffFindPivot(tmp,&f, space->noc) != MTX_NVAL) 
	    break;
    }

    /* Clean up and return the result
       ------------------------------ */
    sysFree(tmp);
    return (i <= 0);
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
