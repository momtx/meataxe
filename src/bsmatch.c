/* ============================= C MeatAxe ==================================
   File:        $Id: bsmatch.c,v 1.1.1.1 2007/09/02 11:06:16 mringe Exp $
   Comment:     Count matching bits.
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */

#include "meataxe.h"

   
/* --------------------------------------------------------------------------
   Local data
   -------------------------------------------------------------------------- */

MTX_DEFINE_FILE_INFO


static int BitCount[256] =  /* Number of '1' bits in binary representation */
{
    0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
};


/**
 ** @addtogroup bs
 ** @{
 **/

/**
 ** Intersection count.
 ** This function calculates the cardinality of the intersection of two bit strings, i.e.,
 ** the number of bits that are set in both @em a and @em b. The arguments must be bit
 ** strings of the same size.
 ** @return Number of bits in the intersection of @em a and @em b, or -1 on error.
 **/

int BsIntersectionCount(const BitString_t *a, const BitString_t *b)
{
    register int i;
    register const unsigned long *ap, *bp;
    int count = 0;

    /* Check the arguments
       ------------------- */
    if (!BsIsValid(a) || !BsIsValid(b))
    {
	MTX_ERROR1("%E",MTX_ERR_BADARG);
	return -1;
    }


    ap = (const unsigned long *) a->Data;
    bp = (const unsigned long *) b->Data;
    for (i = a->BufSize; i > 0; --i)
    {
	register unsigned long a = *ap++ & *bp++;
	while (a != 0)
	{
	    count += BitCount[a & 0xFF];
	    a >>= 8;
	}
    }
    return count;
}
    
  
/**
 ** @}
 **/
