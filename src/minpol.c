/* ============================= C MeatAxe ==================================
   File:        $Id: minpol.c,v 1.2 2007-11-16 00:13:44 mringe Exp $
   Comment:     Minimal polynomial of a matrix
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include <string.h>

MTX_DEFINE_FILE_INFO


/* ------------------------------------------------------------------
   Function prototypes
   ------------------------------------------------------------------ */


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static long fl, nor;		/* Field and size */
static int *APiv = NULL;	/* Pivot table (A) */
static int *CPiv = NULL;	/* Pivot table (C) */
static char *CIsPiv = NULL;	/* Pivot flags */
static PTR mat = NULL,		/* The matrix */
    A = NULL,			/* Work space (spin-up) */
    B = NULL,			/* Work space II (coefficients) */
    C = NULL,			/* Work space III (basis) */
    seed = NULL;
static long CDim;		/* Dimension reached so far */
static long ADim;
static Poly_t *mpol = NULL;

/* ------------------------------------------------------------------
   cleanup() - Free memory
   ------------------------------------------------------------------ */

static void cleanup()

{
    if (mat != NULL) SysFree(mat);
    if (A != NULL) SysFree(A);
    if (B != NULL) SysFree(B);
    if (C != NULL) SysFree(C);
    if (seed != NULL) SysFree(seed);
    if (APiv != NULL) SysFree(APiv);
    if (CPiv != NULL) SysFree(CPiv);
    if (CIsPiv != NULL) SysFree(CIsPiv);
    if (mpol != NULL) PolFree(mpol);
    mpol = NULL;
    mat = A = B = C = NULL;
    CPiv = APiv = NULL;
    CIsPiv = NULL;
}


/* ------------------------------------------------------------------
   init()
   ------------------------------------------------------------------ */

static int init(Matrix_t *matrix)

{	
    /* Set some global variables
       ------------------------- */
    if (matrix->Nor != matrix->Noc) 
	return MTX_ERROR1("%E",MTX_ERR_NOTSQUARE), -1;
    fl = matrix->Field;
    nor = matrix->Nor;
    CDim = 0;

    /* Allocate memory
       --------------- */
    cleanup();
    FfSetField(fl);
    FfSetNoc(nor);
    if ((mat = FfAlloc(nor)) == NULL
	|| (A = FfAlloc(nor+1)) == NULL		/* We need 1 extra row !! */
	|| (B = FfAlloc(nor+1)) == NULL
	|| (C = FfAlloc(nor+1)) == NULL
	|| (seed = FfAlloc(1)) == NULL
	|| (mpol = PolAlloc(fl,0)) == NULL
	|| (APiv = NALLOC(int,nor+2)) == NULL
	|| (CPiv = NALLOC(int,nor+2)) == NULL
	|| (CIsPiv = NALLOC(char,nor+2)) == NULL
       )
	return -1;

    /* Initialize memory
       ----------------- */
    memcpy(mat,matrix->Data,FfCurrentRowSize * nor);
    memset(CIsPiv,0,(size_t)(nor+2));

    return 0;
}


/* ------------------------------------------------------------------
   mkpoly() - Make polynomial for the latest cyclic subspace
   ------------------------------------------------------------------ */

static Poly_t *mkpoly()
{
    int k;
    PTR x;
    Poly_t *pol;
    
    pol = PolAlloc(fl,ADim);
    if (pol == NULL)
	return NULL;
    x = FfGetPtr(B,ADim);
    for (k = 0; k < ADim; ++k)
	pol->Data[k] = FfExtract(x,k);
    pol->Data[ADim] = FF_ONE;
    return pol;
}


/* ------------------------------------------------------------------
   spinup_cyclic() - Spin up one cyclic subaspace
   ------------------------------------------------------------------ */

static void spinup_cyclic()

{   PTR a, b, c;
    long k, pv;
    FEL f;

    a = A;
    b = B;
    c = FfGetPtr(C,CDim);
    FfCopyRow(a,seed);
    FfMulRow(b,FF_ZERO);
    FfInsert(b,0,FF_ONE);
    ADim = 0;
    while ((pv = FfFindPivot(a,&f)) >= 0)
    {	
	PTR x, y;

    	/* Add new vector to basis A
	   ------------------------- */
	FfCopyRow(c,a);
	APiv[ADim++] = pv;
	FfStepPtr(&a);
	FfStepPtr(&b);

    	/* Add new vector to basis C
	   ------------------------- */
	FfCleanRow(c,C,CDim,CPiv);
	if ((pv = FfFindPivot(c,&f)) >= 0)
	{
	    CPiv[CDim++] = pv;
	    CIsPiv[pv] = 1;
	    FfStepPtr(&c);
	}

	/* Calculate the next vector
	   ------------------------- */
	FfMapRow(seed,mat,nor,a);
	FfCopyRow(seed,a);
	FfMulRow(b,FF_ZERO);
	FfInsert(b,ADim,FF_ONE);

	/* Clean with existing basis vectors
	   --------------------------------- */
	x = A; y = B;
	for (k = 0; k < ADim; ++k)
	{   
	    f = FfNeg(FfDiv(FfExtract(a,APiv[k]),FfExtract(x,APiv[k])));
	    FfAddMulRow(a,x,f);
	    FfAddMulRow(b,y,f);
	    FfStepPtr(&x);
	    FfStepPtr(&y);
	}
    }
}


/** @addtogroup charpol
 ** @{
 **/

/**
 ** Minimal Polynomial.
 ** This function returns one factor of the minimal polynomial of
 ** a given matrix. Further calls with a 0 argument return
 ** more factors or 0, if there are no more factors. 
 ** Note that the factors obtained in this way are in general not irreducible.
 **
 ** If @a mat is different from 0, %MinPolFactor() initializes its 
 ** internal data and starts computing one cyclic subspace. Then, the polynomial 
 ** of the matrix restricted to that cyclic subspace is constructed and returned 
 ** to the caller.
 **
 ** If @a mat is 0 on the next call, %MinPolFactor() resumes at 
 ** the point where it returned the last time, calculates the next cyclic 
 ** subspace and so on, until the complete space is exhausted.
 **
 ** @attention Since the function stores information across multiple calls in
 ** static buffers, your program must not use
 ** %MinPolFactor() on more than one matrix at the same time.
 ** @param mat Pointer to the matrix, or 0.
 ** @return One factor of the minmal polynomial or 0 if there are no more factors.
 **/

Poly_t *MinPolFactor(Matrix_t *mat)
{
    int i;
    Poly_t *ggt, *l, *h;

    /* If called with mat != NULL, initialize everything
       ------------------------------------------------- */
    if (mat != NULL)
    {
	if (!MatIsValid(mat))
	    return NULL;
    	if (init(mat) != 0)
	{
    	    MTX_ERROR("Cannot initialize");
	    return NULL;
	}
    }

    /* If there is nothing left to do, return NULL
       ------------------------------------------- */
    /* we could call cleanup() here... */
    if (CDim >= nor) 
	return NULL;

    /* Prepare the next seed vector
       ---------------------------- */
    FfSetField(fl);
    FfSetNoc(nor);
    for (i = 0; i < nor && CIsPiv[i] != 0; ++i);
    MTX_VERIFY(i < nor);
    FfMulRow(seed,FF_ZERO);
    FfInsert(seed,i,FF_ONE);

    /* Spin up and return the polynomial 
       --------------------------------- */
    spinup_cyclic();
    h = mkpoly();    
    ggt = PolGcd(h,mpol);
    l = PolDivMod(h,ggt);
    PolFree(h);
    PolFree(ggt);
    PolMul(mpol,l);
    return l;
}


/**
 ** Minimal polynomial.
 ** This function calculates the minimal polynomial of a matrix in
 ** factored form. The return value is a pointer to a FPoly_t structure
 ** containing the irreducible factors of the minimal polynomial.
 ** @param mat Pointer to the matrix.
 ** @return Minimal polynomial of @a mat, or 0 on error.
 **/

FPoly_t *MinPol(Matrix_t *mat)

{
    Poly_t *f;
    FPoly_t *mp;

    /* Check the argument
       ------------------ */
    if (!MatIsValid(mat))
	return NULL;

    f = MinPolFactor(mat);
    mp = FpAlloc();
    while (f != NULL)
    {
    	FPoly_t *ff = Factorization(f);
    	FpMul(mp,ff);
    	PolFree(f);
    	FpFree(ff);
	f = MinPolFactor(NULL);
    }
    return mp;
}


/**
 ** @}
 **/
