////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Convert between MeatAxe and GAP format.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include <meataxe.h>


/// @addtogroup ff

/// Convert to GAP format.
/// This function takes a field element and returns the GAP
/// representation of this element. The return value is a pointer
/// to a static buffer which is overwritten on each call.
/// @param f Field element.
/// @return Pointer to the GAP representtion of @a f.

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

/// @}
