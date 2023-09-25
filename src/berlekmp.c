////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Berlekamp factorization
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"

/// @addtogroup poly
/// @{


/// Internal representation of a factor.
typedef struct { Poly_t *p; long n; } factor_t;

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Squarefree factorization.
///
/// Returns a factorization of the given polynomial into squarefree factors.
/// The returned list is terminated by an entry with p=NULL, n=0.

static factor_t *factorsquarefree(const Poly_t *pol)
{
    long int j,k,tdeg;
    Poly_t *t, *w, *v;
    FEL *tbuf, *t0buf;

    ffSetField(pol->field);
    uint32_t exp = 0;       // Degree of F over its prime field.
    for (uint32_t ltmp = 1; ltmp != ffChar; ++exp, ltmp *= ffChar);

    Poly_t* t0 = polDup(pol);
    uint32_t e = 1;

    // Allocate the list of factors. There can be at most <pol->Degree>
    // factors, but we need one extra entry to terminate the list.
    size_t nfactors = 0;
    factor_t *factors = NALLOC(factor_t,pol->degree + 1);

    /* Main loop
       --------- */
    while (t0->degree > 0) {
	Poly_t *der = polDerive(polDup(t0));
        t = polGcd(t0,der);
	polFree(der);
        v = polDivMod(t0,t);
        polFree(t0); 
        for (k = 0; v->degree > 0; )
	{
	    Poly_t *tmp;
	    if ( ++k % ffChar == 0 )
	    {
		tmp = polDivMod(t,v);
              	polFree(t);
              	t = tmp;
	      	k++;
	    }
	    w = polGcd(t,v);
	    factors[nfactors].p = polDivMod(v,w);
	    factors[nfactors].n = e * k;
	    if (factors[nfactors].p->degree > 0)
	       	++nfactors;			/* add to output */
	    else
	    	polFree(factors[nfactors].p);	/* discard const. */
            polFree(v);
	    v = w;
	    tmp = polDivMod(t,v);
            polFree(t); 
            t = tmp; 
	} 
	polFree(v);

	/* shrink the polynomial */
      	tdeg = t->degree;
      	e *= ffChar;
      	if ( tdeg % ffChar != 0 )
	    mtxAbort(MTX_HERE,"error in t, degree not divisible by prime");
      	t0 = polAlloc(ffOrder,tdeg/ffChar);
      	t0buf = t0->data;
      	tbuf = t->data;
      	for (j = t0->degree; j >= 0; --j)
	{
	    FEL el, el1;
	    long lexp;

	    el = *tbuf;
	    tbuf += ffChar;
	    el1 = el;
	    /* this is a bug, it should apply the e-1 power of the
	        frobenius K. Lux 16th 11 1996 */

	    /* replace el1 by el1^(ffChar^(exp-1)) */
	    for ( lexp = 0; lexp < exp-1; lexp++ )
	    {
		long n;
		/* replace el1 by el1^ffChar */
		el = el1;
		for ( n = ffChar-1;  0 < n;  n-- ) {
		    el1 = ffMul(el1,el);
		}
	    }
	    *t0buf = el1;
	    ++t0buf;
	}
        polFree(t); 
    }
    polFree(t0);

    // Terminate the list
    MTX_ASSERT(nfactors <= pol->degree);
    factors = NREALLOC(factors, factor_t, nfactors + 1);
    factors[nfactors].p = NULL;
    factors[nfactors].n = 0;
    return factors;
}
		      


/* ------------------------------------------------------------------
   makekernel() - Determines the matrix of the frobenius and returns
   its nullspace.
   ------------------------------------------------------------------ */

static Matrix_t *makekernel(const Poly_t *pol)
{
    FEL* const pbuf = pol->data;
    const long pdeg = pol->degree;
    int k, xshift;
    const long fl = pol->field;
    Matrix_t* materg = matAlloc(fl,pdeg,pdeg);
    PTR rowptr = materg->data;
    MESSAGE(3,("makekernel: fl=%ld pdeg=%ld\n", fl, pdeg));

    FEL* xbuf = NALLOC(FEL,pdeg+1);
    for (k = 0; k <= pdeg; ++k) 
	xbuf[k] = FF_ZERO;
    xbuf[0] = FF_ONE;

    for (k = 0; k < pdeg; ++k)
    {
	int l;
	for (l = 0; l < pdeg; ++l) 
	    ffInsert(rowptr,l,xbuf[l]);
	ffInsert(rowptr,k,ffSub(xbuf[k],FF_ONE));
	ffStepPtr(&rowptr, pdeg);
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
	    f = ffNeg(ffDiv(xbuf[pdeg],pbuf[pdeg]));
	    for (l = pdeg-1; l >= 0; --l)
		xbuf[l] = ffAdd(xbuf[l],ffMul(pbuf[l],f));
	    xbuf[pdeg] = FF_ZERO;
	}
    }
    sysFree(xbuf);
    return matNullSpace__(materg);
 } 


/* ------------------------------------------------------------------
   berlekamp() - Find the irreducible factors of a squarefree polynomial.
   ------------------------------------------------------------------ */

static Poly_t **berlekamp(const Poly_t *pol, const Matrix_t  *kernel)
{
    PTR vec = kernel->data;
    int i;

    Poly_t** list = NALLOC(Poly_t *,kernel->nor+1);
    Poly_t** list2 = NALLOC(Poly_t *,kernel->nor+1);
    list[0] = polDup(pol);
    int nfactors = 1;
    Poly_t* t = polAlloc(kernel->field,kernel->noc-1);

    /* Loop through all kernel vectors */
    for (int j = 2; j <= kernel->nor; ++j)
    {
	int ngcd = 0;
	ffStepPtr(&vec, kernel->noc);			/* Next vector */
	if (nfactors == kernel->nor) 
	    break;	/* Done? */
	for (i = 0; i < kernel->noc; ++i)
	    t->data[i] = ffExtract(vec,i);
	for (i = kernel->noc-1; t->data[i] == FF_ZERO; --i);
	t->degree = i;
	for (i = 0; i < nfactors; )
	{
	    long s;
	    if (list[i]->degree <= 1) {++i; continue;}
	    for (s = 0; s < ffOrder; ++s)
	    {
		Poly_t *gcd;

		t->data[0] = ffFromInt(s);
		MTX_ASSERT(i >= 0 && i < nfactors);
		gcd = polGcd(list[i],t);
		if (gcd->degree >= 1)
		{
		    MTX_ASSERT(ngcd >= 0 && ngcd < kernel->nor+1);
		    list2[ngcd++] = gcd;
		}
		else
		    polFree(gcd);
		if (ngcd == kernel->nor) break;	/* Done? */
	    }
	    if (ngcd > 0)
	    {
		int p;
		MTX_ASSERT(i >= 0 && i < nfactors);
		polFree(list[i]);
		for (p = i; p < nfactors - 1; ++p)
		    list[p] = list[p+1];
		--nfactors;
	    }
	    else
		++i;
	    if (nfactors == kernel->nor) break;		/* Done? */
	}
	if (ngcd > 0)
	{
	    int p;
	    for (p = 0; p < ngcd; ++p)
	    {
		MTX_ASSERT(nfactors >= 0 && nfactors < kernel->nor+1);
	        list[nfactors++] = list2[p];
	    }
	}
    }
    polFree(t);
    sysFree(list2);

    MTX_ASSERT(nfactors >= 0 && nfactors < kernel->nor+1);
    list[nfactors] = NULL;	/* Terminate the list */
    MTX_ASSERT(nfactors == kernel->nor);
    return list;
}




   
/// Factor a polynomial.
/// This function decomposes a polynomial into irreducible factors using the
/// Berlekamp algorithm.
/// @param pol Polynomial to factor.
/// @return The factorization of @a pol.

FPoly_t *Factorization(const Poly_t *pol)
{
    FPoly_t *factors = fpAlloc();

    /* Step 1: Squarefree factorization
       -------------------------------- */
    factor_t *list = factorsquarefree(pol);
    if (list == NULL)
    {
	mtxAbort(MTX_HERE,"Squarefree factorization failed");
    }

    /* Step 2: Decompose the squarefree factors using Berlekamp's algorithm
       -------------------------------------------------------------------- */
    for (factor_t* l = list; l->p != NULL; ++l)
    {
	Matrix_t *kernel;
	Poly_t **irr, **i;

	kernel = makekernel(l->p);
	if ((irr = berlekamp(l->p,kernel)) == NULL)
	{
	    mtxAbort(MTX_HERE,"Berlekamp factorization failed");
	}
	matFree(kernel);
	for (i = irr; *i != NULL; ++i)
	{
	    fpMulP(factors,*i,l->n);
	    polFree(*i);
	}

	/* Clean up
	   -------- */
	sysFree(irr);
	polFree(l->p);
    }

    /* Clean up
       -------- */
    sysFree(list);
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

    p = polAlloc(5,20);
    while (1)
    {
	if (++count % 100 == 0)
	{
	    printf("%d\n",count);
	    fflush(stdout);
	}
	memset(p->data,0,21);
	p->data[20] = FF_ONE;
	for (i = mtxRandomInt(20); i > 0; --i)
	    p->data[RandInt(20)] = ffFromInt(mtxRandomInt(5));
	/*polprint("p",p);*/
	ff = f = polFactorization(p);
	q = polAlloc(5,0);
	for (; f->n != 0; ++f)
	{
    	    int i;
    	    /*polprint("factor",f->p);*/
    	    for (i = 0; i < f->n; ++i) 
		polMul(q,f->p);
    	    /*printf("exp = %ld\n",f->n);*/
	    polFree(f->p);
	}
	free(ff);
	/*polprint("q",q);*/
	if (polCompare(p,q)) break;
	polFree(q);
    }
    polprint("p",p);
    polprint("q",q);
    return 0;
}

#endif
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
