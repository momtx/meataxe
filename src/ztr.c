////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Transpose a matrix.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"



/* ------------------------------------------------------------------
   Global Data
   ------------------------------------------------------------------ */


static const char *iname, *oname;
static int fl;
static PTR m1, m2;
static int nor, noc, j;
static FILE *f;


static MtxApplicationInfo_t AppInfo = { 
"ztr", "Transpose",
"SYNTAX\n"
"    ztr [-QV] <Mat> <Result>"
"\n"
"ARGUMENTS\n"
"    <Mat> ................... Input file name\n"
"    <Result> ................ Output file name\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"\n"
"FILES\n"
"    <Mat> ................... I The matrix\n"
"    <Result> ................ O The transposed matrix\n"
};

static MtxApplication_t *App = NULL;



static int ReadMatrix()

{

    /* Read matrix
       ----------- */
    if ((f = ffReadHeader(iname,&fl,&nor,&noc)) == NULL)
	return 1;
    if (fl < 2) 
    {
	mtxAbort(MTX_HERE,"%s: %s",iname,MTX_ERR_NOTMATRIX);
	return 1;
    }
    ffSetField(fl);
    m1 = ffAlloc(nor, noc);
    ffReadRows(f,m1,nor, noc);
    fclose(f);
    return 0;
}



static int Init(int argc, char **argv)

{
    App = appAlloc(&AppInfo,argc,argv);
    if (App == NULL)
	return -1;

    /* Command line.
       ------------- */
    if (appGetArguments(App,2,2) < 0)
	return -1;
    iname = App->ArgV[0];
    oname = App->ArgV[1];

    if (ReadMatrix() != 0)
	return 1;
    return 0;
}


static void Cleanup()

{
    appFree(App);
}






/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char **argv)

{
    if (Init(argc,argv) != 0)
    {
	mtxAbort(MTX_HERE,"Initialization failed");
	return 1;
    }

    /* Transpose
       --------- */
    m2 = ffAlloc(1, nor);
    if ((f = ffWriteHeader(oname, fl, noc, nor)) == NULL)
    {
	mtxAbort(MTX_HERE,"Cannot open output file");
	return 1;
    }
    for (j = 0; j < noc; ++j)
    {
	ffMulRow(m2, FF_ZERO, nor);	/* Fill with zeros */
	ffExtractColumn(m1, nor, noc, j, m2);
	ffWriteRows(f, m2, 1, nor);
    }
    fclose(f);

    Cleanup();
    return EXIT_OK;
}




/**
@page prog_ztr ztr - Transpose

@section ztr_syntax Command Line
<pre>
ztr [@em Options] @em Mat @em Result
</pre>

@par @em Options
Standard options, see @ref prog_stdopts.

@par @em Mat
    Input matrix.

@par @em Result
    Transposed matrix.

@section ztr_inp Input Files

@par @em Mat
    Input matrix.

@section ztr_out Output Files

@par @em Result
    Transposed matrix.


@section ztr_desc Description
This program transposes a matrix.

*/

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
