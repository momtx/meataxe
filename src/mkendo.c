/* ============================= C MeatAxe ==================================
   File:        $Id: mkendo.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Find endomorphisms
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */

#include "meataxe.h"

MTX_DEFINE_FILE_INFO

/// @defgroup endo Endomorphisms
/// @{

static Matrix_t *MakeEndo(const MatRep_t *rep, const Matrix_t *sb1, 
    const Matrix_t *vec)

{
    Matrix_t *sb2;
    Matrix_t *endo;

    /* Make standard basis from <vec>
       ------------------------------ */
    sb2 = SpinUp(vec,rep,SF_FIRST|SF_CYCLIC|SF_STD,NULL,NULL);
    MTX_VERIFY(sb2 != NULL && sb2->Nor == sb2->Noc);

    /* The linear mapping that maps <sb1> on <Matrix_t> is the endomorphism
       we are looking for!
       --------------------------------------------------------------- */
    endo = MatInverse(sb1);
    MatMul(endo,sb2);
    MatFree(sb2);
    return endo;
}




/// Calculate the endomorphism ring.
/// This function calculates a basis for the endomorphism ring E of an
/// irreducible module. The dimension of E is equal to the degree [F':F]
/// of the splitting field extension for the module.
///
/// The function takes three input arguments: @a ngen is the number of generators
/// of the algebra, and @a gen points to an array of @a ngen matrices containing
/// the representations of the generators. @a nsp is the kernel of an 
/// identifying word for the module, i.e., an algebra element a with
/// dim(ker(a_V))=[F':F].
///
/// On success, basis of E (a set of matrices) is stored into @a endo.
/// If the function fails, no matrices are stored in @a endo. 
/// @param rep
///   Pointer to a matrix representation.
/// @param nsp
///   Idword kernel (see below).
/// @param endo
///   Pointer to result buffer.
/// @return
///   0 on success, -1 on error.

int MakeEndomorphisms(const MatRep_t *rep, const Matrix_t *nsp,	
	Matrix_t *endo[])

{
    Matrix_t *sb1;		/* Standard bases */
    int nendo;			/* # of endomorphisms obtained so far */

    MTX_ASSERT(nsp->Nor > 0);
    MTX_ASSERT(rep->NGen > 0);

    /* Take the first vector from <nsp> and make the standard basis.
       ------------------------------------------------------------- */
    sb1 = SpinUp(nsp,rep,SF_FIRST|SF_CYCLIC|SF_STD,NULL,NULL);
    MTX_VERIFY(sb1 != NULL && sb1->Nor == sb1->Noc);

    /* Take the identity as the first basis element for E
       -------------------------------------------------- */
    endo[0] = MatId(rep->Gen[0]->Field,rep->Gen[0]->Nor);
    nendo = 1;

    /* For each of the remaining vectors v_2,..v_d in <nsp>, construct the
       endomorphism that maps v_1 to v_j.
       -------------------------------------------------------------------- */
    while (nendo < nsp->Nor)
    {
	Matrix_t *vec = MatCutRows(nsp,nendo,1);
	endo[nendo] = MakeEndo(rep,sb1,vec);
	MatFree(vec);
	if (endo[nendo] == NULL)	/* Error */
	    break;
	++nendo;
    }

    /* Clean up
       -------- */
    if (nendo < nsp->Nor)
    {
	while (nendo > 0)
	    MatFree(endo[--nendo]);
    }
    MatFree(sb1); 

    /* Return error code
       ----------------- */
    return nendo == nsp->Nor ? 0 : -1;
}

/// @}
