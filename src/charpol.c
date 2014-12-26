/* ============================= C MeatAxe ==================================
   File:        $Id: charpol.c,v 1.2 2007-11-16 00:13:44 mringe Exp $
   Comment:     Characteristic polynomial of a matrix
   --------------------------------------------------------------------------
   (C) Copyright 1999 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include <string.h>

MTX_DEFINE_FILE_INFO


/* --------------------------------------------------------------------------
   Global data
   -------------------------------------------------------------------------- */

/// @defgroup charpol Characteristic and Minimal Polynomials
/// @{


/// Seed for Characteristic Polynomial.
/// This variable is used by CharPolFactor() to select the first
/// seed vector. By default, CharPolSeed has the value 0, i.e., the
/// first seed vector is (1,0,...,0). Assigning the value 1 selects
/// the start vector (0,1,...,0) in all subsequent calls to
/// CharPolFactor().
/// If CharPolSeed is out of bounds, CharPolFactor() will reset it to 0.

long CharPolSeed = 0;		/* Seed */


/* --------------------------------------------------------------------------
   Local data
   -------------------------------------------------------------------------- */

static long fl, nor;		/* Field and size */
static long *piv = NULL;	/* Pivot table */
static char *ispiv = NULL;	/* Pivot flags */
static PTR mat = NULL,		/* The matrix */
    A = NULL,			/* Work space (for spin-up) */
    B = NULL;			/* Work space II (coefficients) */
static long dim;		/* Dimension reached so far */
static long n;			/* Dimension of cyclic subspace */



/* ------------------------------------------------------------------
   cleanup() - Free memory
   ------------------------------------------------------------------ */

static void cleanup()

{
    if (mat != NULL) SysFree(mat);
    if (A != NULL) SysFree(A);
    if (B != NULL) SysFree(B);
    if (piv != NULL) SysFree(piv);
    if (ispiv != NULL) SysFree(ispiv);
    mat = A = B = NULL;
    piv = NULL;
    ispiv = NULL;
}


/* ------------------------------------------------------------------
   init()
   ------------------------------------------------------------------ */

static int init(const Matrix_t *matrix)

{	
    /* Set some global variables
       ------------------------- */
    if (matrix->Nor != matrix->Noc) 
    {
	MTX_ERROR1("%E",MTX_ERR_NOTSQUARE);
	return -1;
    }
    fl = matrix->Field;
    nor = matrix->Nor;
    dim = 0;
    if (CharPolSeed < 0 || CharPolSeed >= nor)
	CharPolSeed = 0;

    /* Allocate memory
       --------------- */
    cleanup();
    FfSetField(fl);
    FfSetNoc(nor);
    if ((mat = FfAlloc(nor)) == NULL
	|| (A = FfAlloc(nor+1)) == NULL
	|| (B = FfAlloc(nor+1)) == NULL
	|| (piv = NALLOC(long,nor+2)) == NULL
	|| (ispiv = NALLOC(char,nor+2)) == NULL
       )
	return -1;

    /* Initialize memory
       ----------------- */
    memcpy(mat,matrix->Data,FfCurrentRowSize*nor);
    memset(ispiv,0,(size_t)(nor+2));

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
    
    pol = PolAlloc(fl,n);
    x = FfGetPtr(B,n);
    for (k = 0; k < n; ++k)
	pol->Data[k] = FfExtract(x,k);
    pol->Data[n] = FF_ONE;
    return pol;
}


/* --------------------------------------------------------------------------
   spinup_cyclic() - Spin up one cyclic subspace
   -------------------------------------------------------------------------- */

static void spinup_cyclic()

{   PTR a, b;
    long pv, k;
    FEL f;

    a = FfGetPtr(A,dim);
    b = B;
    FfMulRow(b,FF_ZERO);
    n = 0;
    while ((pv = FfFindPivot(a,&f)) >= 0)
    {	
	PTR x, y;

    	/* Add new vector to basis
	   ----------------------- */
	piv[dim+n] = pv;
	ispiv[pv] = 1;
	FfInsert(b,n,FF_ONE);
	++n;

	/* Calculate the next vector
	   ------------------------- */
	x = a;
	FfStepPtr(&a);
	FfMapRow(x,mat,nor,a);
	y = b;
	FfStepPtr(&b);
	FfMulRow(b,FF_ZERO);
	for (k = 0; k < nor-1; ++k)
	    FfInsert(b,k+1,FfExtract(y,k));
    
	/* Clean with existing basis vectors
	   --------------------------------- */
	x = A;
	y = B;
	for (k = 0; k < dim+n; ++k)
	{   
	    f = FfDiv(FfExtract(a,piv[k]),FfExtract(x,piv[k]));
	    FfAddMulRow(a,x,FfNeg(f));
	    if (k >= dim)
	    {	
		FfAddMulRow(b,y,FfNeg(f));
		FfStepPtr(&y);
	    }
	    FfStepPtr(&x);
	}
    }
    dim += n;
}




/// Characteristic Polynomial.
/// This function returns one factor of the characteristic polynomial of
/// a given matrix. Further calls with a 0 argument return
/// more factors or 0, if there are no more factors. 
/// Note that the factors obtained in this way are in general not irreducible.
///
/// Here is how %CharPolFactor() works: If @a mat is different from 0,
/// %CharPolFactor() initializes its internal data and starts
/// computing one cyclic subspace. The choice of starting vector for this
/// first subspace depends on the global variable CharPolSeed.
/// Usually, this variable has a value of 0, corresponding to the vector
/// (1,0,...,0). Then, the polynomial of the matrix restricted to 
/// that cyclic subspace is constructed and returned to the caller.
///
/// If @a mat is 0 on the next call, %CharPolFactor() resumes at 
/// the point where it returned the last time, calculates the next cyclic 
/// subspace and so on, until the complete space is exhausted.

/// @attention Since the function uses static variables to store
/// information across multiple calls, your program must not use
/// %CharPolFactor() on more than one matrix at the same time.
/// @param mat Pointer to the matrix.
/// @return A factor of the characteristic polynomial, or 0 if there are no more factors.

Poly_t *CharPolFactor(const Matrix_t *mat)
{
    PTR seed;
    int i;

    /* If called with mat != NULL, initialize everything
       ------------------------------------------------- */
    if (mat != NULL)
    {
	if (!MatIsValid(mat))
	    return NULL;
    	if (init(mat) != 0)
	{
    	    MTX_ERROR("Initialization failed");
	    return NULL;
	}
    }

    /* If there is nothing left to do, return NULL
       ------------------------------------------- */
    /* we could call cleanup() here... */
    if (dim >= nor) 
	return NULL;

    /* Prepare the next seed vector
       ---------------------------- */
    FfSetField(fl);
    FfSetNoc(nor);
/*    seed = FfGetPtr(A,dim,FfNoc);*/
    seed = (PTR)((char *)A + dim * FfCurrentRowSize);
    if (dim == 0)
	i = CharPolSeed;
    else
    	for (i = 0; i < nor && ispiv[i] != 0; ++i);
    FfMulRow(seed,FF_ZERO);
    FfInsert(seed,i,FF_ONE);

    /* Spin up and return the polynomial 
       --------------------------------- */
    spinup_cyclic();
    return mkpoly();    
}




/// Characteristic Polynomial.
/// This function calculates the characteristic polynomial of a matrix in
/// factored form. The return value is a pointer to a FPoly_t structure
/// containing the irreducible factors of the characteristic polynomial.
/// @param mat Pointer to the matrix.
/// @return The characteristic polynomial of @a mat or 0 on error.

FPoly_t *CharPol(const Matrix_t *mat)
{
    FPoly_t *cpol = FpAlloc();
    Poly_t *p;

    if (cpol == NULL)
    {
	MTX_ERROR("Cannot allocate result");
	return NULL;
    }

    for (p = CharPolFactor(mat); p != NULL; p = CharPolFactor(NULL))
    {
	FPoly_t *fac = Factorization(p);
	if (fac == NULL)
	{
	    MTX_ERROR("Factorization failed");
	    return NULL;
	}
    	PolFree(p);
	FpMul(cpol,fac);
	FpFree(fac);
    }
    return cpol;
}


/// @}
