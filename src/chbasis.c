/* ============================= C MeatAxe ==================================
   File:        $Id: chbasis.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Change basis.
   --------------------------------------------------------------------------
   (C) Copyright 1997 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"

MTX_DEFINE_FILE_INFO


/**
!section obj.matrep
 ** Change basis.
 ** @param rep
    Matrix representation.
 ** @param trans   
    Basis transformation matrix.
 ** @return
    $0$ on success, $-1$ on error.
!description
    This function performs a change of basis on a matrix representation.
    |trans| is the transformation matrix $T$. The rows of |trans| must 
    contain the new basis vectors, expressed in the old basis. Then, the 
    transformed generators are given by
    \[
	g_i' = T g_i T^{-1}
    \]
 ** @see 
 **/

int MrChangeBasis(MatRep_t *rep, const Matrix_t *trans)

{
    Matrix_t *bi;
    int i;

    /* Check arguments
       --------------- */
    if (!MrIsValid(rep))
    {
	MTX_ERROR1("rep: %E",MTX_ERR_BADARG);
	return -1;
    }
    if (!MatIsValid(trans))
    {
	MTX_ERROR1("trans: %E",MTX_ERR_BADARG);
	return -1;
    }
    if (rep->NGen <= 0)
	return 0;
    if (trans->Field != rep->Gen[0]->Field || 
	trans->Nor != rep->Gen[0]->Nor ||
	trans->Noc != rep->Gen[0]->Noc)
    {
	MTX_ERROR1("%E",MTX_ERR_INCOMPAT);
	return -1;
    }


    /* Basis transformation
       -------------------- */
    if ((bi = MatInverse(trans)) == NULL) 
    {
	MTX_ERROR("Basis transformation is singular");
	return -1;
    }
    for (i = 0; i < rep->NGen; ++i)
    {
	Matrix_t *tmp = MatDup(trans);
	MatMul(tmp,rep->Gen[i]);
	MatMul(tmp,bi);
        MatFree(rep->Gen[i]);
	rep->Gen[i] = tmp;
    }
    MatFree(bi);
    return 0;
}



int ChangeBasisOLD(const Matrix_t *M, int ngen, const Matrix_t *gen[],
	Matrix_t *newgen[])

{
    Matrix_t *bi, *tmp;
    int i;

    MTX_VERIFY(ngen >= 0);
    if (!MatIsValid(M))
	return -1;
    if ((bi = MatInverse(M)) == NULL) 
    {
	MTX_ERROR("Matrix is singular");
	return -1;
    }
    for (i = 0; i < ngen; ++i)
    {
	tmp = MatDup(M);
	MatMul(tmp,gen[i]);
	MatMul(tmp,bi);
	if ((const Matrix_t **)newgen == gen)
	    MatFree(newgen[i]);
	newgen[i] = tmp;
    }
    MatFree(bi);
    return 0;
}


