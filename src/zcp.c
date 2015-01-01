////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Characteristic polynomial of a matrix
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"
#include <string.h>
#include <stdlib.h>




/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

MTX_DEFINE_FILE_INFO

static const char *fname;		/* Matrix file name */
static Matrix_t *mat;			/* The matrix */
static int opt_G = 0;			/* GAP output */
static int opt_f = 0;			/* Factorization */
static int opt_m = 0;			/* Minimal polynomial */
static FPoly_t *cpol;			/* Char. polynomial (with -f) */


static MtxApplicationInfo_t AppInfo = { 
"zcp", "Characteristic and Minimal Polynomial",
"SYNTAX\n"
"    zcp [-GQVfm] <File>\n"
"\n"
"ARGUMENTS\n"
"    <File> .................. The matrix\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -G ...................... GAP output\n"
"    -m ...................... Calculate the minimal polynomial\n"
"    -f ...................... Factor the polynomial\n"
"\n"
"FILES\n"
"    <File> .................. I A square matrix\n"
};

static MtxApplication_t *App = NULL;




////////////////////////////////////////////////////////////////////////////////////////////////////

static int Init(int argc, const char **argv)

{	
    if ((App = AppAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;

    /* Parse command line
       ------------------ */
    opt_G = AppGetOption(App,"-G --gap");
    opt_f = AppGetOption(App,"-f --factorize");
    opt_m = AppGetOption(App,"-m --minimal-polynomial");
    if (opt_G) MtxMessageLevel = -100;
    if (AppGetArguments(App,1,1) != 1)
	return -1;
    fname = App->ArgV[0];

    /* Read the matrix
       --------------- */
    mat = MatLoad(fname);
    if (mat == NULL)
	return -1;

    return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static int first;

static void write_init()

{
    first = 1;
    if (opt_m)
    {
    	if (opt_G)
	    printf("MeatAxe.MinPol:=[\n");
    	else
	    printf("MINIMAL POLYNOMIAL:\n");
    }
    else
    {
    	if (opt_G)
	    printf("MeatAxe.CharPol:=[\n");
    	else
	    printf("CHARACTERISTIC POLYNOMIAL:\n");
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void write_end()
{
    int i;

    if (opt_G)
	printf("];\n");
    else
    {
	if (opt_f)
	{
	    for (i = 0; i < cpol->NFactors; ++i)
    	    {
		printf("( ");
		PolPrint(NULL,cpol->Factor[i]);
		printf(" )^%d\n",cpol->Mult[i]);
    	    }
	}
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void write_one(const Poly_t *pol)
{ 
    if (opt_G)
    {
	int i;
	if (!first)
	    printf(",\n");
	printf("[");
	for (i = 0; i < pol->Degree; ++i)
	    printf("%s,",FfToGap(pol->Data[i]));
	printf("%s]",FfToGap(pol->Data[i]));
    }
    else if (!opt_f)
    {
	PolPrint(NULL,pol);
	printf("\n");
    }
    else
    {
	FPoly_t *f = Factorization(pol);
	FpMul(cpol,f);
	FpFree(f);
    }
    first = 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void Cleanup()

{
    if (mat != NULL)
	MatFree(mat);
    if (cpol != NULL)
	FpFree(cpol);
    if (App != NULL)
	AppFree(App);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char **argv)

{	
    Poly_t *p;
    int rc = 0;

    if (Init(argc,argv) != 0)
    {
	MTX_ERROR("Initialization failed");
	return -1;
    }
    if (opt_f) 
	cpol = FpAlloc();
    write_init();
    if (opt_m)
    {
    	for (p = MinPolFactor(mat); p != NULL; p = MinPolFactor(NULL))
	{
	    write_one(p);
	    PolFree(p);
	}
    }
    else
    {
    	for (p = CharPolFactor(mat); p != NULL; p = CharPolFactor(NULL))
	{
	    write_one(p);
	    PolFree(p);
	}
    }
    write_end();


    Cleanup();
    return rc;
}


/**
@page prog_zcp zcp - Characteristic Polynomial

@section zcp_syntax Command Line
<pre>
zcp @em Options [-Gfm] @em Matrix
</pre>

@par @em Options
Standard options, see @ref prog_stdopts

@par -G
Generate GAP output.

@par -f
Factor the polynomial.

@par -m
Calculate the minimal polynomial

@par @em Mat
Input matrix

@section zcp_inp Input Files
@par @em Mat
Input matrix

@section zcp_desc Description

This program reads in a square matrix and calculates its characteristic or
minimal polynomial. With no options, the characteristic polynomial is
computed in a partially factored form (see below). With "-m" the polynomial
is split into irreducible factors. 
Without "-G", the output is in text format. Each line contains one
factor of the characteristic or minimal polynomial.
The "-G" option may be used to generate output which
is readable by the GAP computer program. The output, then, is a
sequence of sequences of finite field elements, representing the
coefficients of the factors in ascending order.

@section zcp_impl Implementation Details
The characteristic polynomial of a matrix A is computed by constructing a sequence 
@f[
        0=U_0<U_1<\ldots<U_n=V
@f]
of A-invariant subspaces, where @f$U_i/U_{i-1}@f$ is A-cyclic.
In the i-th step, @f$U_i@f$ is constructed by chosing a random vector
@f$u\in V\setminus U_{i-1}@f$ and calculating @f$u,uA,uA^2,\ldots@f$ until some
linear combination of these vectors falls into @f$U_{i-1}@f$. This yields
a polynomial @f$p_i(x)@f$ with @f$up_i(A)\in U_{i-1}@f$. @f$p_i(x)@f$ is the
characteristic polynomial of @f$A@f$ on @f$U_i/U_{i-1}@f$, and the full
characteristic polynomial of @f$A@f$ is the product of all @f$p_i@f$'s.

The algorithm for the minimal polynomial uses the same technique of
constructing a sequence @f$(U_i)@f$ of invariant subspaces. In the
i-th step, images @f$uA,uA^2,\ldots@f$ of a seed vector @f$u@f$ are
calculated, until a linear combination of these vectors vanishes
(this is the main difference to the algorithm above). This yields a
polynomial @f$p_i(x)@f$ of minimal degree with @f$up_i(A)=0@f$, and the
minimal polynomial of A is the least common multiple of the
@f$p_i@f$'s.

*/
