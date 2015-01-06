/* ============================= C MeatAxe ==================================
   File:        $Id: c-charpol.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Check characteristic polynomial.
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "c-charpol.h"
#include "meataxe.h"
#include "check.h"

#include <stdlib.h>
#include <stdio.h>

/*MTX_DEFINE_FILE_INFO*/



/* --------------------------------------------------------------------------
   TestCharPol() - Characteristic polynomial
   -------------------------------------------------------------------------- */

static void CheckPoly(Poly_t *p, int degree, ...)

{
    int i;
    va_list al;
    va_start(al,degree);
    if (p->Degree != degree) Error("Bad degree");
    for (i = 0; i <= degree; ++i)
    {
	if (p->Data[i] != FfFromInt(va_arg(al,int)))
	    Error("Bad polynomial");
    }
    va_end(al);
    PolFree(p);
}


void TestCharPol(unsigned flags)
{
    Matrix_t *a;
    Poly_t *pf;

    SelectField(2);
    a = MkMat(6,6,
	1,0,0,0,0,0,  	0,1,1,0,0,0, 	0,0,0,0,1,0, 
	0,0,1,1,0,0,	0,0,0,0,0,1,	0,0,0,0,1,1);

    pf = CharPolFactor(a); CheckPoly(pf,1,1,1);
    pf = CharPolFactor(NULL); CheckPoly(pf,4,0,1,0,0,1);
    pf = CharPolFactor(NULL); CheckPoly(pf,1,1,1);
    pf = CharPolFactor(NULL);
    if (pf != NULL) Error("too many factors");
    MatFree(a);
    flags = 0;
}

