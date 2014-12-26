/* ============================= C MeatAxe ==================================
   File:        $Id: permcmp.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Compare permutations.
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include <string.h>

MTX_DEFINE_FILE_INFO


/// @addtogroup perm
/// @{


/// Compare two permutations.
/// This function compares two permutations. If the permutations are equal, 
/// the return value is 0. Otherwise the return value is positive if @em a 
/// is "greater" than  @em b, and negative if @em a is "less" than @em b. The 
/// ordering of permutations is defined as follows. If the permutations have
/// different degrees, the permutations with the smaller degree is smaller.
/// Otherwise, the result of the comparison is unspecified.
///
/// Note that the return value -1 does not necessarily mean that an error 
/// occured.
/// @param a Pointer to the first permutation.
/// @param b Pointer to the second permutation.
/// @return 0, if the matrices are equal, a nonzero value otherwise, -1 on error.

int PermCompare(const Perm_t *a, const Perm_t *b)
{
    int i;

    /* Check the arguments
       ------------------ */
    if (!PermIsValid(a) || !PermIsValid(b))
    {
	MTX_ERROR1("%E",MTX_ERR_BADARG);
	return -1;
    }

    /* Compare degrees
       --------------- */
    if ((i = a->Degree - b->Degree) != 0)
	return i;

    /* Compare the entries
       ------------------- */
    i = memcmp(a->Data,b->Data,sizeof(long) * a->Degree);
    if (i < 0)
	return -1;
    else if (i > 0)
	return 1;
    return 0;
}

/// @}
