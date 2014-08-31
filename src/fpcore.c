/* ============================= C MeatAxe ==================================
   File:        $Id: fpcore.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Basic factored polynomial functions.
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

#define FP_MAGIC 0x17B69244


/** 
 ** @addtogroup poly
 ** @{
 **/

/** @class FPoly_t
 ** A Factored Polynomial.
 ** This structure contains a polynomial which is split into factors. The factors
 ** need not be irreducible.
 **/
   

/**
 ** Check a Factored Polynomial.
 ** @param p The polynomial.
 ** @return 1 if @em p is a valid factores polynomial, 0 otherwise.
 **/

int FpIsValid(const FPoly_t *p)
{
    int i;
    if (p == NULL)
    {
	MTX_ERROR("NULL polynomial");
	return 0;
    }
    if (p->Magic != FP_MAGIC || p->NFactors < 0 || p->BufSize < p->NFactors)
    {
	MTX_ERROR3("Invalid FPoly_t: Magic=%d, NFactors=%d, MaxLen=%d",
	    (int)p->Magic,p->NFactors,p->BufSize);
	return 0;
    }
    if (p->Factor == NULL || p->Mult == NULL)	
    {
	MTX_ERROR2("Invalid FPoly_t: Factor:%s, Mult:%s",
	    p->Factor == 0 ? "NULL":"ok",
	    p->Mult == 0 ? "NULL":"ok");
	return 0;
    }
    for (i = 0; i < p->NFactors; ++i)
    {
	if (!PolIsValid(p->Factor[i]))
	{
	    MTX_ERROR("Invalid factor");
	    return 0;
	}
	if (p->Mult[i] < 0)
	{
	    MTX_ERROR1("Invalid multiplicity %d",p->Mult[i]);
	    return 0;
	}
	if (i > 0 && p->Factor[i]->Field != p->Factor[0]->Field)
	{
	    MTX_ERROR("Factors over different fields");
	    return 0;
	}
    }
    return 1;
}



/**
 ** Allocate a Factored Polynomial.
 ** This function creates a new Fpoly_t structure.
 ** The new polynomial is empty, i.e., it has no factors.
 ** @return Pointer to the new FPoly_t structure or 0 on error.
 **/

FPoly_t *FpAlloc()

{
    FPoly_t *x;

    x = ALLOC(FPoly_t);
    if (x == NULL)
    {
	MTX_ERROR1("%E",MTX_ERR_NOMEM);
	return NULL;
    }
    x->BufSize = 5;
    x->Factor = NALLOC(Poly_t *,x->BufSize);
    if (x->Factor == NULL)
    {
	SysFree(x);
	MTX_ERROR1("%E",MTX_ERR_NOMEM);
	return NULL;
    }
    x->Mult = NALLOC(int,x->BufSize);
    if (x->Mult == NULL)
    {
	SysFree(x->Factor);
	SysFree(x);
	MTX_ERROR1("%E",MTX_ERR_NOMEM);
	return NULL;
    }
    x->NFactors = 0;
    x->Magic = FP_MAGIC;
    return x;
}




/**
 ** Free a Factored Polynomial.
 ** @return 0 on success, -1 on error.
 ** @see FPoly_t FpAlloc
 **/

int FpFree(FPoly_t *x)
{
    int i;

    /* Check the argument
       ------------------ */
    if (!FpIsValid(x))
	return -1;

    /* Free all factors
       ---------------- */
    for (i = 0; i < x->NFactors; ++i)
	PolFree(x->Factor[i]);

    /* Free the <FPoly_t> structure
       ---------------------------- */
    SysFree(x->Factor);
    SysFree(x->Mult);
    memset(x,0,sizeof(FPoly_t));
    SysFree(x);
    return 0;
}


/** 
 ** @}
 **/
