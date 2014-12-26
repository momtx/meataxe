/* ============================= C MeatAxe ==================================
   File:        $Id: poladd.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Add two polynomials.
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



/* ------------------------------------------------------------------
   resize() - Extend buffer.
   ------------------------------------------------------------------ */

static int resize(Poly_t *p, int newdeg)
{
    int i;
    FEL *x;

    if (p->Degree < newdeg)
    {
	if (p->BufSize < newdeg + 1)	/* Allocate new buffer */
    	{
	    x = NREALLOC(p->Data,FEL,newdeg + 1);
	    if (x == NULL)
	    {
		MTX_ERROR("Cannot extend polynomial");
		return 1;
	    }
	    p->Data = x;
	    p->BufSize = newdeg + 1;
    	}
    	for (i = p->Degree + 1; i <= newdeg; ++i)
	    p->Data[i] = FF_ZERO;
        p->Degree = newdeg;
    }
    return 0;
}





/// Add Polynomials.
/// Thus function adds @em src to @em dest.
/// The polynomials must be over the same field.
/// @param dest Pointer to the first polynomial.
/// @param src Pointer to the second polynomial.
/// @return @em dest, or 0 on error.

Poly_t *PolAdd(Poly_t *dest, const Poly_t *src)
{
    FEL *s, *d;
    int i;

    if (!PolIsValid(src) || !PolIsValid(dest))
	return NULL;
    if (dest->Field != src->Field) 
    {
	MTX_ERROR1("%E",MTX_ERR_INCOMPAT);
	return NULL;
    }
    if ((i = src->Degree) == -1) 
	return dest;	/* src = 0 */

    FfSetField(src->Field);
    if (resize(dest,i)) 
    {
	MTX_ERROR("Cannot resize: %S");
	return NULL;
    }
    s = src->Data;
    d = dest->Data;
    for (; i >= 0; --i)
    {
       *d = FfAdd(*d,*s++);
       ++d;
    }
    Pol_Normalize(dest);
    return dest;
}

/// @}
