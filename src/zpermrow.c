/* ============================= C MeatAxe ==================================
   File:        $Id: zpermrow.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Map a vector under a permutation.
   --------------------------------------------------------------------------
   (C) Copyright 1997 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"

#ifdef PARANOID
MTX_DEFINE_FILE_INFO
#endif



/**
 ** @addtogroup ff2
 ** @{
 **/
 
/**
 ** Multiply a Vector by a Permutation.
 ** This function multiplies the vector @a row from the right with the permutation @a perm
 ** and stores the result in @a result. Multiplication of vectors by permutations is defined
 ** as follows: if the permutation maps i to k, then the i-ith mark of the vector is stored
 ** in the k-th position of the result. 
 **
 ** Note: @a result and @a row must not overlap. Otherwise the result is undefined.
 **/

void FfPermRow(PTR row, const long *perm, PTR result)
{
    register FEL f;
    register int i;
    register const long *p = (long *)perm;

#ifdef PARANOID
    if (row == result)
	MTX_ERROR("row = result: undefined result!");
#endif

    for (i = 0; i < FfNoc; ++i)
    {
#ifdef PARANOID
	if (*p < 0 || *p > FfNoc)
	    MTX_ERROR2("Invalid point %d in permutation, noc=%d",(int)*p,FfNoc);
#endif
	f = FfExtract(row,i);
	FfInsert(result,*p++,f);
    }
}

/**
 ** @}
 **/
