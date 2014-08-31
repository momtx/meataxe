/* ============================= C MeatAxe ==================================
   File:        $Id: bsdup.c,v 1.1.1.1 2007/09/02 11:06:16 mringe Exp $
   Comment:     Duplicate a bit string.
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
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
 ** @addtogroup bs
 ** @{
 **/

/**
 ** Duplicate a bit string.
 ** Returns an independent copy of the argument, or 0 on error.
 **/

BitString_t *BsDup(const BitString_t *src)
{
    BitString_t *n;

    if (!BsIsValid(src))
    {
	MTX_ERROR1("%E",MTX_ERR_BADARG);
	return NULL;
    }
    n = BsAlloc(src->Size);
    if (n == NULL)
	return NULL;
    memcpy(n->Data,src->Data,src->BufSize * sizeof(long));
    return n;
}

/**
 ** @}
 **/
