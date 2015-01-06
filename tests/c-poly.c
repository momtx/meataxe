/* ============================= C MeatAxe ==================================
   File:        $Id: c-poly.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Check functions for matrices.
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "c-poly.h"
#include "meataxe.h"
#include "check.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

MTX_DEFINE_FILE_INFO

static int ErrorFlag = 0;

static void MyErrorHandler(const MtxErrorRecord_t *err)
{
    ErrorFlag = 1;
    err = NULL;
}

static int CheckError()
{
    int i = ErrorFlag;
    ErrorFlag = 0;
    return i;
}


Poly_t *RndPol(int fl, int mindeg, int maxdeg)
{
    Poly_t *p;
    int deg;
    int i;

    deg = MtxRandomInt(maxdeg-mindeg+1) + mindeg;
    p = PolAlloc(fl,deg);
    for (i = 0; i <= deg-1; ++i)
	p->Data[i] = FfFromInt(MtxRandomInt(FfOrder));
    return p;
}


/* --------------------------------------------------------------------------
   TestPolAlloc() - Permutation allocation
   -------------------------------------------------------------------------- */

#define NPOLY 5

void TestPolAlloc1(int fl)

{
    static const int deg[NPOLY] = { -1,0,5,10,200 };
    Poly_t *p[NPOLY];
    MtxErrorHandler_t *old_err_handler;
    int i;

    for (i = 0; i < NPOLY; ++i) p[i] = PolAlloc(fl,deg[i]);
    for (i = 0; i < NPOLY; ++i) 
    {
	int k;
	PolIsValid(p[i]);
	MTX_VERIFY(p[i]->Degree == deg[i]);
	for (k = 0; k < p[i]->Degree; ++k)
	    MTX_VERIFY(p[i]->Data[k] == FF_ZERO);
	if (p[i]->Degree >= 0)
	    MTX_VERIFY(p[i]->Data[p[i]->Degree] == FF_ONE);
    }
    for (i = 0; i < NPOLY; ++i) 
	if (PolFree(p[i]) != 0) Error("PolFree() failed");
    old_err_handler = MtxSetErrorHandler(MyErrorHandler);
    for (i = 0; i < NPOLY; ++i) 
	if (PolIsValid(p[i]) || !CheckError()) Error("PolIsValid() failed");
    MtxSetErrorHandler(old_err_handler);
}




void TestPolAlloc(unsigned flags)
{
    while (NextField() > 0)
	TestPolAlloc1(FfOrder);
    flags = 0;
}



/* --------------------------------------------------------------------------
   TestPolAlloc() - Permutation allocation
   -------------------------------------------------------------------------- */

#define NPOLY 5

static void TestPolCompare2(int fl1, int deg1, int fl2, int deg2, int result)
{
    Poly_t *a, *b;
    int i;

    a = PolAlloc(fl1,deg1);
    b = PolAlloc(fl2,deg2);
    i = PolCompare(a,b);
    if (i != result)
	Error("PolCompare(a,b) failed: result=%d, expected %d",i,result);
    i = PolCompare(b,a);
    if (i != -result)
	Error("PolCompare(b,a) differs from -PolCompare(a,b)");
    PolFree(a);
    PolFree(b);
}



static void TestPolCompare1(int fl)

{
    Poly_t *a, *b;
    const int deg = 10;
    int i;

    a = PolAlloc(fl,deg);
    b = PolAlloc(fl,deg);
    if (PolCompare(a,b) != 0)
	Error("PolCmp(a,a) != 0");
    for (i = 0; i < deg; ++i)
    {
	a->Data[i] = FF_ONE;
	if (PolCompare(a,b) == 0)
	    Error("PolCmp() == 0 on different polnomials");
	b->Data[i] = FF_ONE;
        if (PolCompare(a,b) != 0)
	    Error("PolCmp(a,a) != 0");
    }
    if (FfGen != FF_ONE)
    {
	a->Data[deg] = FfGen;
	if (PolCompare(a,b) == 0)
	    Error("PolCmp() == 0 on different polnomials");
    }
    PolFree(a);
    PolFree(b);
}




void TestPolCompare(unsigned flags)
{
    TestPolCompare2(2,-1,3,-1,-1);
    TestPolCompare2(2, 0,3,-1,-1);
    TestPolCompare2(2,10,3, 0,-1);
    TestPolCompare2(2,0,2,-1,1);
    TestPolCompare2(2,10,2,9,1);
    TestPolCompare2(3,0,3,0,0);
    TestPolCompare2(3,-1,3,-1,0);

    while (NextField() > 0)
	TestPolCompare1(FfOrder);
    flags = 0;
}



/* --------------------------------------------------------------------------
   TestPolAdd() - Permutation addition
   -------------------------------------------------------------------------- */

static void TestPolAdd1(int fl)

{
    Poly_t *a, *b, *c;

    a = PolAlloc(fl,-1);
    b = PolAlloc(fl,0);
    c = PolAlloc(fl,0);
    PolAdd(a,b);
    if (PolCompare(a,c) != 0) Error("PolAdd(a,0)!=a failed");
    PolFree(a);
    PolFree(b);
    PolFree(c);
}

static void TestPolAdd2(int fl)

{
    Poly_t *a, *b;
    int i;

    a = PolAlloc(fl,-1);
    for (i = -1; i <= 10; ++i)
    {
	b = PolAlloc(fl,i);
	PolAdd(a,b);
	PolFree(b);
    }
    if (a->Degree != 10)
	Error("Wrong degree");
    for (i = 0; i <= 10; ++i)
    {
	if (a->Data[i] != FF_ONE)
	    Error("Wrong coefficient");
    }
    PolFree(a);
}



void TestPolAdd(unsigned flags)
{
    while (NextField() > 0)
    {
	TestPolAdd1(FfOrder);
	TestPolAdd2(FfOrder);
    }
    flags = 0;
}



/* --------------------------------------------------------------------------
   TestPolMul() - Permutation multiplication
   -------------------------------------------------------------------------- */


static void TestPolMul1(int fl)

{
    Poly_t *a, *b;

    /* Check x * 0 = 0 
       --------------- */
    a = PolAlloc(fl,1);
    b = PolAlloc(fl,-1);
    PolMul(a,b);
    if (a->Degree != -1) Error("x * 0 != 0");
    PolFree(a);
    PolFree(b);

    /* Check x * 1 = x 
       --------------- */
    a = PolAlloc(fl,1);
    b = PolAlloc(fl,0);
    PolMul(a,b);
    if (a->Degree != 1 || a->Data[0] != FF_ZERO || a->Data[1] != FF_ONE) 
	Error("x * 1 != x");
    PolFree(a);
    PolFree(b);

    /* Check (x+1)*(x^3-x) = x^4+x^3-x^2-x 
       ----------------------------------- */
    a = PolAlloc(fl,1);
    b = PolAlloc(fl,3);
    a->Data[0] = FF_ONE;
    b->Data[1] = FfNeg(FF_ONE);
    PolMul(a,b);
    if (a->Degree != 4 || a->Data[0] != FF_ZERO || a->Data[1] != FfNeg(FF_ONE)
	|| a->Data[2] != FfNeg(FF_ONE) || a->Data[3] != FF_ONE 
	|| a->Data[4] != FF_ONE ) 
	Error("Error in (x+1)(x^3-2)");
    PolFree(a);
    PolFree(b);
}



static void TestPolMul2(int fl)
{
    int i;
    MtxRandomInit(fl);
    for (i = 0; i < 20; ++i)
    {
	Poly_t *a, *b, *c, *ab, *ba;
	a = RndPol(fl,0,100);
	b = RndPol(fl,0,100);
	c = RndPol(fl,0,100);

	ab = PolDup(a);
	PolMul(ab,b);
	ba = PolDup(b);
	PolMul(ba,a);
	if (PolCompare(ab,ba) != 0) Error("ab != ba");

	PolMul(ab,c);
	PolMul(b,c);
	PolMul(a,b);
	if (PolCompare(a,ab) != 0) Error("(ab)c != a(bc)");

	PolFree(a);
	PolFree(b);
	PolFree(c);
	PolFree(ab);
	PolFree(ba);
    }
}

void TestPolMul(unsigned flags)
{
    while (NextField() > 0)
    {
	TestPolMul1(FfOrder);
	TestPolMul2(FfOrder);
    }
    flags = 0;
}

/* --------------------------------------------------------------------------
   TestPolMul() - Permutation multiplication
   -------------------------------------------------------------------------- */

static void TestPolGcd1(int fl)

{
    int n;

    MtxRandomInit(fl);
    for (n = 100; n > 0; --n)
    {
	Poly_t *a, *b, *gcd, *result[3];

	/* Create two random polynomials
	   ----------------------------- */
	a = RndPol(fl,0,100);
	b = RndPol(fl,0,100);

	/* Calculate the g.c.d.
	   -------------------- */
	if ((gcd = PolGcd(a,b)) == NULL)
	    Error("PolGcd");
	if (PolGcdEx(a,b,result) != 0)
	    Error("PolGcdEx");

	/* Check the result
	   ---------------- */
	if (PolCompare(gcd,result[0]) != 0)
	    Error("PolGcd != PolGcdEx");
	PolMul(result[1],a);
	PolMul(result[2],b);
	PolAdd(result[1],result[2]);
	if (PolCompare(gcd,result[1]) != 0)
	    Error("PolGcdEx coeffs");
	if (PolMod(a,gcd) == NULL)
	    Error("PolMod");
	if (PolMod(b,gcd) == NULL)
	    Error("PolMod");
	if (a->Degree >= 0 || b->Degree >= 0)
	    Error("gcd");
	PolFree(a);
	PolFree(b);
	PolFree(gcd);
	PolFree(result[0]);
	PolFree(result[1]);
	PolFree(result[2]);
    }
}



void TestPolGcd(unsigned flags)
{
    while (NextField() > 0)
    {
	TestPolGcd1(FfOrder);
    }
    flags = 0;
}
