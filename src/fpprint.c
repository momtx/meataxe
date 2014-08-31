/* ============================= C MeatAxe ==================================
   File:        $Id: fpprint.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Print a factored polynomial.
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */

#include "meataxe.h"

   
/* --------------------------------------------------------------------------
   Local data
   -------------------------------------------------------------------------- */

/*MTX_DEFINE_FILE_INFO*/



/** 
 ** @addtogroup poly
 ** @{
 **/

/**
 ** Print a factored polynomial.
 ** This function prints a factored polynomial to the standard output.
 ** If @em name is not 0, "name=" is printed before the polynomial.
 ** @param name Name of the polynomial or 0.
 ** @param p Pointer to the factored polynomial.
 ** @return 0 on success, -1 on error.
 **/

int FpPrint(const char *name, const FPoly_t *p)
{
    int i;
    if (!FpIsValid(p))
	return -1;
    if (name != NULL) 
	printf("%s =",name);
    for (i = 0; i < p->NFactors; ++i)
    {   
	int e = p->Mult[i];
	if (i > 0) printf("    * ");
	printf("(");
	PolPrint(NULL,p->Factor[i]);
	if (e > 1)
	    printf(")^%d\n",e);
	else
	    printf(")\n");
    }
    return 0;
}

/** 
 ** @}
 **/

