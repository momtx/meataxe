////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Split a representation
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>

/// @addtogroup g_spinup
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

static int CheckArguments(const Matrix_t *subspace, const MatRep_t *rep)
{
    mrValidate(MTX_HERE, rep);
    matValidate(MTX_HERE, subspace);
    if (subspace->pivotTable == NULL)
    {
	mtxAbort(MTX_HERE,"%s",MTX_ERR_NOTECH);
	return -1;
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Split a Representation.
/// Given a matrix representation of an algebra A and an A-invariant 
/// subspace U, this function calculates two new matrix representations 
/// corrsponding to the subspace and quotient, respectively.
///
/// @a subspace is a basis for the invariant subspace. This matrix must be in
/// echelon form. @a rep is the representation to be split. @a sub and @a quot
/// are pointers to variables, where the representation on subspace and 
/// quotient, respectively, will be stored. Both @a sub and @a quot can be 0
/// if the corresponding representation is not needed. If @a sub is
/// not 0, <tt>*sub</tt> must be 0 or a pointer to a valid matrix 
/// representation, which will be destroyed before the result is stored 
/// in @a *sub|. The same remark applies to @a quot.
///
/// The function checks that the subspace is indeed invariant under the
/// given representation. However, this check is carried out only when
/// the subspace is calculated, i.e., if @a sub is not 0. The function
/// also check if subspace and representation are compatible. If any of 
/// these checks fails, the return value is -1.
///
/// Internally, Split() uses SAction() and QAction() to calculate
/// the action of the generators on the subspace and quotient.
///
/// The following example shows how to use SpinUp() to find an invariant subspace.
/// If a proper subspace is found, the representation is split.
/// @code
/// MatRep_t *Rep;
/// Matrix_t *seed;
/// Matrix_t *subspace;
/// ...
/// subspace = SpinUp(seed,Rep,SF_FIRST|SF_SUB);
/// if (subspace->nor > 0 && subspace->nor < subspace->noc)
/// {
///     MatRep_t *Sub = NULL, *Quot = NULL;
///     printf("Split!\n");
///     Split(subspace,Rep,&Sub,&Quot);
/// }
/// @endcode
///
/// @see SpinUp QAction SAction
/// @param subspace Pointer to an invariant subspace.
/// @param rep Matrix representation to split.
/// @param sub Matrix representation on the subspace.
/// @param quot Matrix representation on the quotient.
/// @return 0 on success, -1 on error.

int Split(const Matrix_t *subspace, const MatRep_t *rep, MatRep_t **sub, MatRep_t **quot)
{
    int g;

    /* Check Arguments.
       ---------------- */
    if (CheckArguments(subspace,rep) != 0)
    {
	mtxAbort(MTX_HERE,"%s",MTX_ERR_BADARG);
	return -1;
    }

    /* Subspace
       -------- */
    if (sub != NULL)
    {
	if (*sub != NULL)
	    mrFree(*sub);
	*sub = mrAlloc(0,NULL,0);
	if (*sub == NULL)
	    return -1;
	for (g = 0; g < rep->NGen; ++g)
	{
	    Matrix_t *gen = SAction(subspace,rep->Gen[g]);
	    mrAddGenerator(*sub,gen,0);
	}
    }

    /* Quotient
       -------- */
    if (quot != NULL)
    {
	if (*quot != NULL)
	    mrFree(*quot);
	*quot = mrAlloc(0,NULL,0);
	if (*quot == NULL)
	    return -1;
	for (g = 0; g < rep->NGen; ++g)
	{
	    Matrix_t *gen = QAction(subspace,rep->Gen[g]);
	    mrAddGenerator(*quot,gen,0);
	}
    }

    return 0;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
