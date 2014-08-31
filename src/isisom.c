/* ============================= C MeatAxe ==================================
   File:        $Id: isisom.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Isomorphism test.
   --------------------------------------------------------------------------
   (C) Copyright 1999 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"


MTX_DEFINE_FILE_INFO



static int CheckArgs(int ngen, Matrix_t  **gen1, const CfInfo *info1,
    Matrix_t **gen2, int use_pw)

{
    int j;

    /* Check arguments
       --------------- */
    MTX_VERIFY(ngen > 0);
    for (j = 0; j < ngen; ++j)
    {
        if (!MatIsValid(gen1[j]) || !MatIsValid(gen2[j]))
	    return -1;
	if (gen1[j]->Nor != gen1[j]->Noc)
	{
	    MTX_ERROR1("gen1[%d]: Matrix not square",j);
	    return -1;
	}
	if (gen2[j]->Nor != gen2[j]->Noc)
	{
	    MTX_ERROR1("gen2[%d]: Matrix not square",j);
	    return -1;
	}
	if (gen1[j]->Field != gen1[0]->Field || gen1[j]->Nor != gen1[0]->Nor)
	{
	    MTX_ERROR1("gen1[%d]: Incompatible matrix",j);
	    return -1;
	}
	if (gen2[j]->Field != gen1[0]->Field)
	{
	    MTX_ERROR1("gen2[%d]: Incompatible matrix",j);
	    return -1;
	}
    }

    if (info1->dim != gen1[0]->Nor)
    {
	MTX_ERROR("Inconsistent cfinfo data");
	return -1;
    }
    if (use_pw && info1->peakword == 0)
    {
	MTX_ERROR("No peak word available");
	return -1;
    }
    if (!use_pw && info1->idword == 0)
    {
	MTX_ERROR("No id word available");
	return -1;
    }
    return 0;
}



/**
!section algo.etc
 ** Compare representations.
 ** @param rep1
    The first representation.
 ** @param info1
    Pointer to an |CfInfo| structure for the first representation.
 ** @param rep2
    The second representation.
 ** @param trans
    Buffer for basis transformation matrix, or |NULL|.
 ** @param use_pw
    If different from zero, use peak word instead of the identifying word 
    (see below).
 ** @return
    1 if the representations are isomorphic, 0 if they are not isomorphic, and
    $-1$ on error.
!description
    This function decides if two irreducible representations are isomorphic.
    |rep1| and |rep2| must be two matrix representations over the same field,
    and with the same number of generators. Furthermore,
    to compare the representations, the function needs an identifying word
    for the first representation, i.e., the fields |info1->idword|, 
    |info1->idpol| and |info1->spl| must be set, and the generators in 
    |rep1| must be in standard basis with respect to the identifying word.
    If |use_pw| is nonzero, the peak word is used instead of the idword.
    In this case, |rep1| must of course be in standard basis with respect 
    to the peak word.

    If the representations are isomorphic, and |trans| is not |NULL|, the
    basis transformation which makes the second representation identical
    to the first is stored into |*trans|. To be more precise, if $g_i$ is 
    the representation of the i-th generator in representation |rep1|,
    $h_i$ in representation |rep2|, and $T$ the 
    matrix returned in |trans|, then $Th_i T^{-1}=g_i$.
 ** @see 
 **/

int IsIsomorphic(const MatRep_t *rep1, const CfInfo *info1,
    const MatRep_t *rep2, Matrix_t  **trans, int use_pw)

{
    int j;
    WgData_t *wg;
    Matrix_t  *word, *m, *seed, *b, *bi;
    int result;

    if (CheckArgs(rep1->NGen,rep1->Gen,info1,rep2->Gen,use_pw) != 0)
	return -1;
 
    /* Check if the dimensions are equal
       --------------------------------- */
    if (rep1->Gen[0]->Nor != rep2->Gen[0]->Nor)
	return 0;

    /* Make the idword on representation 2
       ----------------------------------- */
    wg = WgAlloc(rep2);
    word = WgMakeWord(wg,use_pw ? info1->peakword : info1->idword);
    m = MatInsert(word,use_pw ? info1->peakpol : info1->idpol);
    MatFree(word);
    WgFree(wg);
    seed = MatNullSpace__(m);
    if (seed->Nor != info1->spl)
    {	
	MatFree(seed);
	return 0;
    }

    /* Make the standard basis
       ----------------------- */
    b = SpinUp(seed,rep2,SF_FIRST|SF_CYCLIC|SF_STD,NULL,NULL);
    MatFree(seed);
    if (b->Nor != b->Noc)
    {
	MatFree(b);
	return 0;
    }
    bi = MatInverse(b);

    /* Compare generators
       ------------------ */
    for (j = 0, result = 0; result == 0 && j < rep2->NGen; ++j)
    {
	Matrix_t *g = MatDup(b);
	MatMul(g,rep2->Gen[j]);
	MatMul(g,bi);
	if (MatCompare(g,rep1->Gen[j]) != 0)
	    result = 1;
	MatFree(g);
    }

    /* Clean up 
       -------- */
    if (trans != NULL && result == 0)
	*trans = b;
    else
	MatFree(b);
    MatFree(bi);

    return (result == 0);
}
