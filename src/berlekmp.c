/* ============================= C MeatAxe ==================================
   File:        $Id: berlekmp.c,v 1.1.1.1 2007/09/02 11:06:16 mringe Exp $
   Comment:     Berlekamp factorization
   --------------------------------------------------------------------------
   (C) Copyright 1997 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"

MTX_DEFINE_FILE_INFO

/// @addtogroup poly
/// @{


/// Internal representation of a factor.
typedef struct { Poly_t *p; long n; } factor_t;

/* ------------------------------------------------------------------
   factorsquarefree() - Squarefree factorization.

   This function factors a polynomial into square-free, but not
   necessarily irreducible, factors. The factors are stored in a
   list of <factor_t> structures. A pointer to this list is passed
   back to the caller. The list is terminated by a ``factor'' with
   p=NULL, n=0.

   Return: List of factors.
   ------------------------------------------------------------------ */

static factor_t *factorsquarefree(const Poly_t *pol)
{
    long int e,j,k,ltmp,tdeg,exp;
    Poly_t *t0, *t, *w, *v;
    FEL *tbuf, *t0buf;
    factor_t *list;
    int nfactors = 0;

    /* Initialize variables
       -------------------- */
    FfSetField(pol->Field);
    for (exp = 0, ltmp = FfOrder; ltmp % FfChar == 0; ++exp, ltmp /= FfChar);
    t0 = PolDup(pol);
    e = 1;

    /* Allocate the list of factors. The can be at most <pol->Degree>
       factors, but we need one extra entry to terminate the list.
       ----------------------------------------------------------- */
    list = NALLOC(factor_t,pol->Degree + 1);

    /* Main loop
       --------- */
    while (t0->Degree > 0)
    {
	Poly_t *der = PolDerive(PolDup(t0));
        t = PolGcd(t0,der);
	PolFree(der);
        v = PolDivMod(t0,t);
        PolFree(t0); 
        for (k = 0; v->Degree > 0; )
	{
	    Poly_t *tmp;
	    if ( ++k % FfChar == 0 )
	    {
		tmp = PolDivMod(t,v);
              	PolFree(t);
              	t = tmp;
	      	k++;
	    }
	    w = PolGcd(t,v);
	    list[nfactors].p = PolDivMod(v,w);
	    list[nfactors].n = e * k;
	    if (list[nfactors].p->Degree > 0)
	       	++nfactors;			/* add to output */
	    else
	    	PolFree(list[nfactors].p);	/* discard const. */
            PolFree(v);
	    v = w;
	    tmp = PolDivMod(t,v);
            PolFree(t); 
            t = tmp; 
	} 
	PolFree(v);

	/* shrink the polynomial */
      	tdeg = t->Degree;
      	e *= FfChar;
      	if ( tdeg % FfChar != 0 )
	    printf("error in t, degree not div. by prime \n");
      	t0 = PolAlloc(FfOrder,tdeg/FfChar);
      	t0buf = t0->Data;
      	tbuf = t->Data;
      	for (j = t0->Degree; j >= 0; --j)
	{
	    FEL el, el1;
	    long lexp;

	    el = *tbuf;
	    tbuf += FfChar;
	    el1 = el;
	    /* this is a bug, it should apply the e-1 power of the
	        frobenius K. Lux 16th 11 1996 */

	    /* replace el1 by el1^(FfChar^(exp-1)) */
	    for ( lexp = 0; lexp < exp-1; lexp++ )
	    {
		long n;
		/* replace el1 by el1^FfChar */
		el = el1;
		for ( n = FfChar-1;  0 < n;  n-- ) {
		    el1 = FfMul(el1,el);
		}
	    }
	    *t0buf = el1;
	    ++t0buf;
	}
        PolFree(t); 
    }
    PolFree(t0);

    /* Terminate the list
       ------------------ */
    list[nfactors].p = NULL;
    list[nfactors].n = 0;
    return list;
}
		      


/* ------------------------------------------------------------------
   makekernel() - Determines the matrix of the frobenius and returns
   its nullspace.
   ------------------------------------------------------------------ */

static Matrix_t *makekernel(const Poly_t *pol)
{
    Matrix_t *materg;
    PTR rowptr;
    FEL *xbuf, *pbuf = pol->Data;
    long pdeg = pol->Degree;
    int k, xshift;
    long fl = pol->Field;

    materg = MatAlloc(fl,pdeg,pdeg);
    rowptr = materg->Data;

    xbuf = NALLOC(FEL,pdeg+1);
    for (k = 0; k <= pdeg; ++k) 
	xbuf[k] = FF_ZERO;
    xbuf[0] = FF_ONE;

    for (k = 0; k < pdeg; ++k)
    {
	int l;
	for (l = 0; l < pdeg; ++l) 
	    FfInsert(rowptr,l,xbuf[l]);
	FfInsert(rowptr,k,FfSub(xbuf[k],FF_ONE));
	FfStepPtr(&rowptr);
        for (xshift = (int) fl; xshift > 0; )
	{
	    FEL f;
	    int d;

	    /* Find leading pos */
	    for (l = pdeg-1; xbuf[l] == FF_ZERO && l >= 0; --l);

	    /* Shift left as much as possible */
	    if ((d = pdeg - l) > xshift) d = xshift;
	    for (; l >= 0; l--) xbuf[l+d] = xbuf[l];
	    for (l = d-1; l >= 0; --l) xbuf[l] = FF_ZERO;
	    xshift -= d;
	    if (xbuf[pdeg] == FF_ZERO) continue;

	    /* Reduce with pol */
	    f = FfNeg(FfDiv(xbuf[pdeg],pbuf[pdeg]));
	    for (l = pdeg-1; l >= 0; --l)
		xbuf[l] = FfAdd(xbuf[l],FfMul(pbuf[l],f));
	    xbuf[pdeg] = FF_ZERO;
	}
    }
    SysFree(xbuf);
    return MatNullSpace__(materg);
 } 


/* ------------------------------------------------------------------
   berlekamp() - Find the irreducible factors of a squarefree polynomial.
   ------------------------------------------------------------------ */

static Poly_t **berlekamp(const Poly_t *pol, const Matrix_t  *kernel)
{
    Poly_t **list, **list2;
    int nfactors;
    PTR vec = kernel->Data;
    int i, j;
    Poly_t *t;

    list = NALLOC(Poly_t *,kernel->Nor+1);
    list2 = NALLOC(Poly_t *,kernel->Nor+1);
    list[0] = PolDup(pol);
    nfactors = 1;
    t = PolAlloc(kernel->Field,kernel->Noc-1);
    FfSetNoc(kernel->Noc);

    /* Loop through all kernel vectors */
    for (j = 2; j <= kernel->Nor; ++j)
    {
	int ngcd = 0;
	FfStepPtr(&vec);			/* Next vector */
	if (nfactors == kernel->Nor) 
	    break;	/* Done? */
	for (i = 0; i < kernel->Noc; ++i)
	    t->Data[i] = FfExtract(vec,i);
	for (i = kernel->Noc-1; t->Data[i] == FF_ZERO; --i);
	t->Degree = i;
	for (i = 0; i < nfactors; )
	{
	    long s;
	    if (list[i]->Degree <= 1) {++i; continue;}
	    for (s = 0; s < FfOrder; ++s)
	    {
		Poly_t *gcd;

		t->Data[0] = FfFromInt(s);
		gcd = PolGcd(list[i],t);
		if (gcd->Degree >= 1)
		{
		    list2[ngcd++] = gcd;
		}
		else
		    PolFree(gcd);
		if (nfactors == kernel->Nor) break;	/* Done? */
	    }
	    if (ngcd > 0)
	    {
		int p;
		PolFree(list[i]);
		for (p = i; p < nfactors -1; ++p)
		    list[p] = list[p+1];
		--nfactors;
	    }
	    else
		++i;
	    if (nfactors == kernel->Nor) break;		/* Done? */
	}
	if (ngcd > 0)
	{
	    int p;
	    for (p = 0; p < ngcd; ++p)
	    {
	        list[nfactors++] = list2[p];
	    }
	}
    }
    PolFree(t);
    SysFree(list2);

    list[nfactors] = NULL;	/* Terminate the list */
    MTX_VERIFY(nfactors == kernel->Nor);
    return list;
}




   
/// Factor a polynomial.
/// This function decomposes a polynomial into irreducible factors using the
/// Berlekamp algorithm.
/// @param pol Polynomial to factor.
/// @return The factorization of @a pol or 0 on error.

FPoly_t *Factorization(const Poly_t *pol)
{
    factor_t *list, *l;
    FPoly_t *factors;    

    /* Allocate result
       --------------- */
    if ((factors = FpAlloc()) == NULL)
    {
	MTX_ERROR("Cannot allocate result");
	return NULL;
    }

    /* Step 1: Squarefree factorization
       -------------------------------- */
    if ((list = factorsquarefree(pol)) == NULL)
    {
	MTX_ERROR("Squarefree factorization failed");
	return NULL;
    }

    /* Step 2: Decompose the squarefree factors using Berlekamp's algorithm
       -------------------------------------------------------------------- */
    for (l = list; l->p != NULL; ++l)
    {
	Matrix_t *kernel;
	Poly_t **irr, **i;

	kernel = makekernel(l->p);
	if ((irr = berlekamp(l->p,kernel)) == NULL)
	{
	    MTX_ERROR("Berlekamp factorization failed");
	    return NULL;
	}
	MatFree(kernel);
	for (i = irr; *i != NULL; ++i)
	{
	    FpMulP(factors,*i,l->n);
	    PolFree(*i);
	}

	/* Clean up
	   -------- */
	SysFree(irr);
	PolFree(l->p);
    }

    /* Clean up
       -------- */
    SysFree(list);
    return factors;
}


/// @}







#if 0

Poly_t *p, *q, *r;
factor_t *f, *ff;

int main(void)
{
    int i;
    int count = 0;

    p = PolAlloc(5,20);
    while (1)
    {
	if (++count % 100 == 0)
	{
	    printf("%d\n",count);
	    fflush(stdout);
	}
	memset(p->Data,0,21);
	p->Data[20] = FF_ONE;
	for (i = MtxRandomInt(20); i > 0; --i)
	    p->Data[RandInt(20)] = FfFromInt(MtxRandomInt(5));
	/*polprint("p",p);*/
	ff = f = PolFactorization(p);
	q = PolAlloc(5,0);
	for (; f->n != 0; ++f)
	{
    	    int i;
    	    /*polprint("factor",f->p);*/
    	    for (i = 0; i < f->n; ++i) 
		PolMul(q,f->p);
    	    /*printf("exp = %ld\n",f->n);*/
	    PolFree(f->p);
	}
	free(ff);
	/*polprint("q",q);*/
	if (PolCompare(p,q)) break;
	PolFree(q);
    }
    polprint("p",p);
    polprint("q",q);
    return 0;
}

#endif
