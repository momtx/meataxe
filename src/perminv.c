/* ============================= C MeatAxe ==================================
   File:        $Id: perminv.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Inverse of a permutation.
   --------------------------------------------------------------------------
   (C) Copyright 1999 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */

#include "meataxe.h"
#include <string.h>

   
/* --------------------------------------------------------------------------
   Local data
   -------------------------------------------------------------------------- */

MTX_DEFINE_FILE_INFO


/**
 ** @addtogroup perm
 ** @{
 **/

/**
 ** Inverse of a permutation
 ** This function calulates the inverse of a permutation.
 ** @param src Pointer to the permutation.
 ** @return The inverse of @em src, or 0 on error.
 **/

Perm_t *PermInverse(const Perm_t *src)
{
    register long i;
    register long *d, *s;
    Perm_t *inv;

    /* Check arguments 
       --------------- */
    if (!PermIsValid(src))
	return NULL;

    /* Allocate a new permutation.
       --------------------------- */
    inv = PermAlloc(src->Degree);
    if (inv == NULL)
    {
	MTX_ERROR("Cannot allocate result buffer");
	return NULL;
    }
    d = inv->Data;
    s = src->Data;

    /* Invert.
       ------- */
    for (i = src->Degree - 1, s = src->Data + i; i >= 0; --i, --s)
	d[*s] = i;

    return inv;
}

/**
 ** @}
 **/
