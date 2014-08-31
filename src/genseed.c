/* ============================= C MeatAxe ==================================
   File:        $Id: genseed.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Seed vector generator.
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */

#include "meataxe.h"
#include <limits.h>
#include <stdlib.h>


/*{{{ -----------------------------------------------------------------------
!definesection algo.pseed
    The seed vector generator is used to walk through the one-dimensional
    subspaces of a given vector space $V$, the {\em seed space}. For each 
    one-dimensional subspace $U<V$ the generator produces a representant
    $u\in U$. These vectors are called {\em seed vectors}.

    Once a basis for the seed space is fixed, each vector of the seed space 
    can be identified by a natural number as follows. Let $b_1$,\ldots,$b_n$ 
    be the basis and $v=\lambda_1b_1+\ldots\lambda_nb_n\in V$ a vector.
    The coefficients $\lambda_i\in\GF{q}$ are identified with
    positive integers in the usual way (see |FfToInt()| for details) and
    \[
	N(v):=\lambda_1+\lambda_2q+\lambda_3q^2+\ldots\lambda_nq^{n-1}
    \]
    is the `number of $v$'. Seed vectors are those vectors where the
    leading coefficient is $1$.
 **/




/* ------------------------------------------------------------------
   Local data
   ------------------------------------------------------------------ */

MTX_DEFINE_FILE_INFO



/**
!section algo.pseed
 ** Make a seed vector.
 ** @param basis
    Seed space basis.
 ** @param lastno
    Previous seed vector (0 to start from beginning).
 ** @param vec
    Buffer for the seed vector.
 ** @return
    Seed vector number, $-1$ on error.
!description
    This function generates the next seed vector. |basis| is a basis of the
    seed space, it need not be in echelon form. |lastno| is the number of the 
    previous seed vector. A value of $0$ starts with the first seed vector.
    The seed vector is stored in |vec|, and its number is returned. If the
    rows of |basis| are not linearly independent, there will be redundant seed
    vectors, but no error occurs.

    Typically, |MakeSeedVector()| is called repeatedly until the whole
    seed space is exhausted. Here is s short example:
    \begin{verbatim}
	long v;
	PTR vec = FfAlloc(1);
	v = MakeSeedVector(basis,0,ptr); 
	while (v > 0)
	{
	    ...
	    v = MakeSeedVector(basis,v,ptr)
	}
    \end{verbatim}
!seealso
 **/

long MakeSeedVector(const Matrix_t *basis, long lastno, PTR vec)

{
    long nextno, x, i;
    int j;

    if (!MatIsValid(basis))
	return -1;
    if (vec == NULL || lastno < 0)
    {
	MTX_ERROR1("%E",MTX_ERR_BADARG);
	return -1;
    }

    /* Find the next seed vector number
       -------------------------------- */
    nextno = lastno + 1;
    for (x = 1; (i = nextno/x) >= FfOrder; x *= FfOrder);
    if (i != 1) 
	nextno = x * FfOrder;

    /* Make the seed vector
       -------------------- */
    FfSetField(basis->Field);
    FfSetNoc(basis->Noc);
    FfMulRow(vec,FF_ZERO);
    for (j = 0, x = nextno; x != 0 && j < basis->Nor; ++j, x /= FfOrder)
    {
	FEL co = FfFromInt(x % FfOrder);
	if (co != FF_ZERO)
	    FfAddMulRow(vec,MatGetPtr(basis,j),co);
    }

    /* Check for overflow
       ------------------ */
    if (x != 0)
	return -1;

    return nextno;
}


