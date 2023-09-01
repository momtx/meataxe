////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Change basis.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"



/// @addtogroup mrep 
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Change basis.
/// This function performs a change of basis on a matrix representation.
/// The transformed generators are given by g'<sub>i</sub> = T g<sub>i</sub> T<sup>-1</sup>.
/// @param rep Matrix representation.
/// @param trans Transformation matrix mapping the old basis to the new basis. In other words,
///    the rows of the matrix are the new basis vectors
/// @return 0 on success, -1 on error.

int mrChangeBasis(MatRep_t *rep, const Matrix_t *trans)
{
    Matrix_t *bi;
    int i;

    /* Check arguments
       --------------- */
    if (!mrIsValid(rep))
    {
	mtxAbort(MTX_HERE,"rep: %s",MTX_ERR_BADARG);
	return -1;
    }
    matValidate(MTX_HERE, trans);

    if (rep->NGen <= 0)
	return 0;
    if (trans->Field != rep->Gen[0]->Field || 
	trans->Nor != rep->Gen[0]->Nor ||
	trans->Noc != rep->Gen[0]->Noc)
    {
	mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
	return -1;
    }


    /* Basis transformation
       -------------------- */
    if ((bi = matInverse(trans)) == NULL) 
    {
	mtxAbort(MTX_HERE,"Basis transformation is singular");
	return -1;
    }
    for (i = 0; i < rep->NGen; ++i)
    {
	Matrix_t *tmp = matDup(trans);
	if (tmp == NULL) {
	    return -1;
	}
	matMul(tmp,rep->Gen[i]);
	matMul(tmp,bi);
        matFree(rep->Gen[i]);
	rep->Gen[i] = tmp;
    }
    matFree(bi);
    return 0;
}



int ChangeBasisOLD(const Matrix_t *M, int ngen, const Matrix_t *gen[],
	Matrix_t *newgen[])
{
    Matrix_t *bi, *tmp;
    int i;

    MTX_ASSERT(ngen >= 0);
    matValidate(MTX_HERE, M);
    if ((bi = matInverse(M)) == NULL) 
    {
	mtxAbort(MTX_HERE,"Matrix is singular");
	return -1;
    }
    for (i = 0; i < ngen; ++i)
    {
	tmp = matDup(M);
	matMul(tmp,gen[i]);
	matMul(tmp,bi);
	if ((const Matrix_t **)newgen == gen)
	    matFree(newgen[i]);
	newgen[i] = tmp;
    }
    matFree(bi);
    return 0;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
