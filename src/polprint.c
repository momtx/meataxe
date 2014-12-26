/* ============================= C MeatAxe ==================================
   File:        $Id: polprint.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Print a polynomial.
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


/// @addtogroup poly
/// @{

/// Print a Polynomial
/// This function prints a polynomial on the standard output in a human-readable
/// form. If @em name is not 0, the name followed by an equal sign is printed 
/// before the polynomial. For example, the statement <tt>PolPrint("P",P)</tt> could
/// produce the following output:
/// <pre>
/// P=3x^2+x+1</pre>
/// @param name Name to print before the polynomial or 0.
/// @param p Pointer to the polynomial.
/// @return 0 on success, -1 on error.

void PolPrint(char *name, const Poly_t *p)
{
    int i,flag = 0;

    if (!PolIsValid(p))
    {
	MTX_ERROR1("%E",MTX_ERR_BADARG);
	return;
    }
    if (name != NULL) 
	printf("%s=",name);
    FfSetField(p->Field);
    if (p->Degree == -1) 
    {
	printf("0x^0");
    }
    for (i = p->Degree; i >= 0; i--)
    {
	if (p->Data[i] != FF_ZERO)
	{
	    if (flag) 
		printf("+");
	    if (p->Data[i] != FF_ONE || i == 0)
		printf("%d",FfToInt(p->Data[i]));
	    if (i > 1) 
		printf("x^%d",i);
	    else if (i == 1) 
		printf("x");
	    flag=1;
       	}
    }
    if (name != NULL) 
	printf("\n");
}




/// @}
