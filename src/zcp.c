////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Characteristic polynomial of a matrix
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"
#include <string.h>
#include <stdlib.h>

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
"    -f ...................... factor the polynomial\n"
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
//    if (opt_G) MtxMessageLevel = -100;
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
	    printf(i == 0 ? "%s" : ",%s",gapFelToString(pol->data[i]));
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
   for (int i = 0; i < fpol->nFactors; ++i) {
      const Poly_t* const factor = fpol->factor[i];
      const int exp = fpol->mult[i];

      if (opt_G) {
         if (!first)
            printf(",\n");
         first = 0;
         printf("[" "["); // Doxygen chokes on two opening brackets
         for (int i = 0; i < factor->degree; ++i)
            printf("%s,",gapFelToString(factor->data[i]));
         printf("%s], %d]",gapFelToString(factor->data[i]), exp);
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
       Charpol_t* state = charpolStart(mat, opt_m ? PM_MINPOL : PM_CHARPOL, 0);
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

////////////////////////////////////////////////////////////////////////////////////////////////////

// *INDENT-OFF*

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
factor the polynomial.

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
0=U₀<U₁< ‥  <Uₙ=V of A-invariant subspaces, where U<sub>i</sub>/U<sub>i-1</sub> is A-cyclic.
In each step, U<sub>i</sub> is constructed by chosing a random vector u∈V∖U<sub>i-1</sub>
and calculating uA, uA², uA³,… until some linear combination of these vectors falls into
U<sub>i-1</sub>. This yields a factor of the characteristic polynomial corresponding to 
U<sub>i</sub>/U<sub>i-1</sub>, and the full characteristic polynomial is the product of
all factors.

The algorithm for the minimal polynomial uses the same technique of constructing a sequence
(U<sub>i</sub>) of invariant subspaces. In each step images uA, uA², uA³,… of a seed vector
u are calculated until some linear combination of the images is zero (this is the main
difference to the algorithm above). This yields a polynomial p<sub>i</sub>(x) of minimal degree
with up<sub>i</sub>(A)=0, and the minimal polynomial of A is the least common multiple of the
p<sub>i</sub>'s.

*/

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
