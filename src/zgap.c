/* ============================= C MeatAxe ==================================
   File:        $Id: zgap.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Convert between MeatAxe and GAP format.
   --------------------------------------------------------------------------
   (C) Copyright 1999 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"


/**
 ** @addtogroup ff2
 **/

/**
 ** Convert to GAP format.
 ** This function takes a field element and returns the GAP
 ** representation of this element. The return value is a pointer
 ** to a static buffer which is overwritten on each call.
 ** @param f Field element.
 ** @return Pointer to the GAP representtion of @a f.
 **/

const char *FfToGap(FEL f)
{
    static char buffer[40];

    if (FfChar == FfOrder)	/* Prime field */
    {
	FEL f2 = FF_ZERO;
   	int k = 0;
    	while (f2 != f)
    	{  
	    f2 = FfAdd(f2,FfGen);
	    ++k;
	}
	sprintf(buffer,"%d*Z(%d)",k,FfOrder);
    }
    else		/* Other field */
    {
	if (f == FF_ZERO)
	    sprintf(buffer,"0*Z(%d)",FfOrder);
	else
	{
	    FEL f2 = FfGen;
	    int k = 1;
	    while (f2 != f)
	    {   
		f2 = FfMul(f2,FfGen);
		++k;
	    }
	    sprintf(buffer,"Z(%d)^%d",FfOrder,k);
	}
    }
    return buffer;
}

/**
 ** @}
 **/
