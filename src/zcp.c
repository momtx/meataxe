////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Characteristic polynomial of a matrix
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"
#include <string.h>
#include <stdlib.h>




/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */


static const char *fname;		/* Matrix file name */
static Matrix_t *mat;			/* The matrix */
static int opt_G = 0;			/* GAP output */
static int opt_f = 0;			/* Factorize the result */
static int opt_m = 0;			/* Minimal polynomial */
static int opt_p = 0;			/* Return a single polynomial */


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
"    -p ...................... Do not factorize, print a single polynomial\n"
"\n"
"FILES\n"
"    <File> .................. I A square matrix\n"
};

static MtxApplication_t *App = NULL;
static int first = 1;



////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char **argv)
{	
    App = appAlloc(&AppInfo,argc,argv);

    // Parse command line
    opt_G = appGetOption(App,"-G --gap");
    opt_f = appGetOption(App,"-f --factorize");
    opt_m = appGetOption(App,"-m --minimal-polynomial");
    opt_p = appGetOption(App,"-p --single-polynomial");
    if (opt_f && opt_p) {
       mtxAbort(MTX_HERE, "-f and -p cannot be combined");
    }
    if (opt_G) MtxMessageLevel = -100;
    appGetArguments(App,1,1);
    fname = App->argV[0];

    // Read the matrix
    mat = matLoad(fname);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void writeBegin()
{
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

static void writeEnd()
{
   if (opt_G)
      printf("];\n");
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void writeP(const Poly_t *pol)
{ 
    if (opt_G) {
	if (!first)
	    printf(",\n");
        first = 0;
	printf("[");
	for (int i = 0; i <= pol->degree; ++i)
	    printf(i == 0 ? "%s" : ",%s",ffToGap(pol->data[i]));
	printf("]");
    }
    else {
	polPrint(NULL,pol);
	printf("\n");
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void writeF(const FPoly_t *fpol)
{ 
   for (int i = 0; i < fpol->NFactors; ++i) {
      const Poly_t* const factor = fpol->Factor[i];
      const int exp = fpol->Mult[i];
      if (opt_G) {
         if (!first)
            printf(",\n");
         first = 0;
         printf("[[");
         for (int i = 0; i < factor->degree; ++i)
            printf("%s,",ffToGap(factor->data[i]));
         printf("%s], %d]",ffToGap(factor->data[i]), exp);
      } else {
         printf("(");
         polPrint(NULL,factor);
         printf(")^%d\n",exp);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void Cleanup()
{
    if (mat != NULL)
	matFree(mat);
    if (App != NULL)
	appFree(App);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{	
    init(argc,argv);
    writeBegin();

    if (!opt_f) {
       // Partial factorization (default) or single polynomial
       Charpol_t* state = charpolAlloc(mat, opt_m ? PM_MINPOL : PM_CHARPOL, 0);
       Poly_t *poly = polAlloc(mat->field, 0);
       Poly_t *factor;
       while ((factor = charpolFactor(state)) != NULL) {
          if (opt_p)
             polMul(poly, factor);
          else
             writeP(factor);
          polFree(factor);
       }
       if (opt_p)
          writeP(poly);
       polFree(poly);
    } else if (opt_f) {
       // Full factorization
       FPoly_t* poly = opt_m ? minpol(mat) : charpol(mat);
       writeF(poly);
    }

    writeEnd();
    Cleanup();
    return 0;
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
minimal polynomial.

With no options, the characteristic polynomial is computed in a partially factored
form (see below). Factors are in general reducible, and the same factor may occur
multiple times at random positions.

With "-f" the polynomial is split into irreducible factors.

With "-p" all factor are combined and a single polynomial is printed.

By default the output is in text format with each factor on a separate line.
The "-G" option may be used to generate output which is readable by the GAP program.
In GAP mode, the is a sequence of sequences of finite field elements, representing the
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
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
