////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Bit Strings, incidence relation
//
// (C) Copyright 1998-2014 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO


/// @addtogroup bs
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Bit string incidence relation.
/// This function returns 1 if and only if every bit which is set in @em a is also set in 
/// @em b. Both bit strings must have the same size.
/// @return 1 if a⊆b, 0 if a⊈b, -1 on error.

int BsIsSub(const BitString_t *a, const BitString_t *b)
{	
    register int i;
    register const unsigned long *ap, *bp;

    /* Check the arguments
       ------------------- */
    if (!BsIsValid(a))
    {
	MTX_ERROR1("a: %E",MTX_ERR_BADARG);
	return -1;
    }
    if (!BsIsValid(b))
    {
	MTX_ERROR1("b: %E",MTX_ERR_BADARG);
	return -1;
    }
    if (a->Size != b->Size)
    {
	MTX_ERROR1("%E",MTX_ERR_INCOMPAT);
	return -1;
    }

    /* Calculate the result
       -------------------- */
    ap = (unsigned long const *) a->Data;
    bp = (unsigned long const *) b->Data;
    for (i = a->BufSize; i > 0; --i, ++ap, ++bp)
    {
	if ((*ap ^ (*ap & *bp)) != 0)
	    return 0;
    }
    return 1;
}

/// @}
