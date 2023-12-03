////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Reduce to (normalized) echelon form
////////////////////////////////////////////////////////////////////////////////////////////////////



#include "meataxe.h"
#include <stdlib.h>




/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */


static int opt_G = 0;	/* GAP output */
static const char *iname, *oname;
static Matrix_t *Mat = NULL;

static MtxApplicationInfo_t AppInfo = { 
"zef", "Echelon Form",
"SYNTAX\n"
"    zef [-GQV] <Inp> <Out>\n"
"\n"
"ARGUMENTS\n"
"    <Inp> ................... Matrix file name\n"
"    <Out> ................... Output file name\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -G ...................... GAP output (implies -Q)\n"
"\n"
"FILES\n"
"    <Inp> ................... I The matrix\n"
"    <Out> ................... I The reduced matrix\n"
};

static MtxApplication_t *App = NULL;





static int Init(int argc, char **argv)

{
    /* Process command line options and arguments.
       ------------------------------------------- */
    if ((App = appAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    opt_G = appGetOption(App,"-G --gap");
//    if (opt_G) MtxMessageLevel = -100;
    if (appGetArguments(App,2,2) != 2)
	return -1;
    iname = App->argV[0];
    oname = App->argV[1];

    return 0;
}



static void Cleanup()

{
    if (App != NULL)
	appFree(App);
    if (Mat != NULL)
	matFree(Mat);
}





/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char **argv)

{
    int rc = 0;

    if (Init(argc,argv) != 0)
    {
	mtxAbort(MTX_HERE,"Initialization failed");
	return 1;
    }

    Mat = matLoad(iname);
    if (Mat == NULL)
	rc = -1;
    if (rc == 0)
    {
	if (matEchelonize(Mat) < 0)
	    rc = 1;
    }
    if (rc == 0)
    {
	/* Normalize.
	   ---------- */
	int i;
	for (i = 0; i < Mat->nor; ++i)
	{
	    PTR rp = matGetPtr(Mat,i);
	    ffMulRow(rp,ffInv(ffExtract(rp,Mat->pivotTable[i])),Mat->noc);
	}
    }
    if (rc == 0)
    {
	matSave(Mat,oname);
    }
    if (rc == 0)
    {
	if (opt_G)
	    printf("MeatAxe.Rank := %d;\n",Mat->nor);
	else
	    MTX_LOGI("RANK %d",Mat->nor);
    }
    Cleanup();
    return rc;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// *INDENT-OFF*

/**
@page prog_zef zef - Echelon Form

@section zef_syntax Command Line
<pre>
zef @em Options [-G] @em Inp @em Out
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts.
@par -G
  GAP output.
@par @em Inp
  Input matrix.
@par @em Out
  Output matrix.

@section zef_inp Input Files
@par @em Inp
  Input matrix.

@section zef_out Output Files
@par @em Out
  Output matrix.

@section zef_desc Description
This program reads in a matrix, performs a Gaussian elimination to put
the matrix into semi-echelon form, and writes out the result.

A matrix is in semi-echelon form, if the first non-zero entry in each
row entry is a 1, and all entries exactly below that 1 are zero. Here
is an example:
<pre>
000112312312342412123
012012312312231231233
000000001345121223233
102012330332333312212
001021230333323123123
000001230311212112121
</pre>
A matrix is in (full) echelon form if the rows are further permuted so
that the first non-zero entry in each row is later than the first
non-zero entry in all previous rows. There is no real need to do this
row permutation in the MeatAxe system, and the permutation is not
done by this program.

At the end, the matrix may have fewer rows than it started with, since
the rows may have been linearly dependent. The rows of the output matrix
are always linearly independent and span the same space
as the rows of the input matrix. A message
<pre>
RANK nnn</pre>
is printed at the end of the run to notify the operator of the size of
the output matrix. This program may be used to find the rank of 
matrix, being faster than the null-space program. There is no need for
the input matrix to be square, and the output may also not be square.

*/
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
