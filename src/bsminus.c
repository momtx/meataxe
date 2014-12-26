////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Bit Strings, difference
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
/// Difference of two bit strings.
/// This function computes the (set theoretical) difference of two bit strings, i.e., a bit
/// in the result is set if and only if it is set in @em dest but not in @em src.
/// The result is stored in @em dest and overwrites the previous value.
/// The arguments must be bit strings of the same size.
/// @return 0 on success, -1 on error.

int BsMinus(BitString_t *dest, const BitString_t *src)
{
    register int i;
    register unsigned long *dp;
    register const unsigned long *sp;

    /* Check the arguments
       ------------------- */
    if (!BsIsValid(dest))
    {
	MTX_ERROR1("dest: %E",MTX_ERR_BADARG);
	return -1;
    }
    if (!BsIsValid(src))
    {
	MTX_ERROR1("src: %E",MTX_ERR_BADARG);
	return -1;
    }
    if (dest->Size != src->Size)
    {
	MTX_ERROR1("%E",MTX_ERR_INCOMPAT);
	return -1;
    }

    /* AND operation
       ------------- */
    dp = (unsigned long *) dest->Data;
    sp = (unsigned long const *) src->Data;
    for (i = src->BufSize; i > 0; --i)
	*dp++ &= ~*sp++;

    return 0;
}


/// @}
