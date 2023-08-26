////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Convert between MeatAxe and GAP format.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"


/// @addtogroup ff
/// @{

/// Convert to GAP format.
/// This function takes a field element and returns the GAP
/// representation of this element. The return value is a pointer
/// to a static buffer which is overwritten on each call.
/// @param f Field element.
/// @return Pointer to the GAP representtion of @a f.

const char *ffToGap(FEL f)
{
    static char buffer[40];

    if (ffChar == ffOrder)	/* Prime field */
    {
	FEL f2 = FF_ZERO;
   	int k = 0;
    	while (f2 != f)
    	{  
	    f2 = ffAdd(f2,ffGen);
	    ++k;
	}
	sprintf(buffer,"%d*Z(%d)",k,ffOrder);
    }
    else		/* Other field */
    {
	if (f == FF_ZERO)
	    sprintf(buffer,"0*Z(%d)",ffOrder);
	else
	{
	    FEL f2 = ffGen;
	    int k = 1;
	    while (f2 != f)
	    {   
		f2 = ffMul(f2,ffGen);
		++k;
	    }
	    sprintf(buffer,"Z(%d)^%d",ffOrder,k);
	}
    }
    return buffer;
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
