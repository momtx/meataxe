////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Transpose a matrix.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"



/* ------------------------------------------------------------------
   Global Data
   ------------------------------------------------------------------ */

MTX_DEFINE_FILE_INFO

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
    if ((f = FfReadHeader(iname,&fl,&nor,&noc)) == NULL)
	return 1;
    if (fl < 2) 
    {
	MTX_ERROR2("%s: %E",iname,MTX_ERR_NOTMATRIX);
	return 1;
    }
    FfSetField(fl);
    FfSetNoc(noc);
    m1 = FfAlloc(nor);
    FfReadRows(f,m1,nor);
    fclose(f);
    return 0;
}



static int Init(int argc, const char **argv)

{
    App = AppAlloc(&AppInfo,argc,argv);
    if (App == NULL)
	return -1;

    /* Command line.
       ------------- */
    if (AppGetArguments(App,2,2) < 0)
	return -1;
    iname = App->ArgV[0];
    oname = App->ArgV[1];

    if (ReadMatrix() != 0)
	return 1;
    return 0;
}


static void Cleanup()

{
    AppFree(App);
}






/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, const char **argv)

{
    if (Init(argc,argv) != 0)
    {
	MTX_ERROR("Initialization failed");
	return 1;
    }

    /* Transpose
       --------- */
    FfSetNoc(nor);
    m2 = FfAlloc(1);
    if ((f = FfWriteHeader(oname,fl,noc,nor)) == NULL)
    {
	MTX_ERROR("Cannot open output file");
	return 1;
    }
    for (j = 0; j < noc; ++j)
    {
	FfSetNoc(nor);
	FfMulRow(m2,FF_ZERO);	/* Fill with zeros */
	FfSetNoc(noc);
	FfExtractColumn(m1,nor,j,m2);
	FfSetNoc(nor);
	FfWriteRows(f,m2,1);
    }
    fclose(f);

    Cleanup();
    return (EXIT_OK);
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

