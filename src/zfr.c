////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Apply the Frobenius automorphism to a matrix.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include <meataxe.h>


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

MTX_DEFINE_FILE_INFO

static MtxApplicationInfo_t AppInfo = { 
"zfr", "Frobenius Automorphism",
"SYNTAX\n"
"    zfr [-QV] <Matrix> <Result>\n"
"\n"
"ARGUMENTS\n"
"    <Matrix> ................ Input file name\n"
"    <Result> ................ Output file name\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"\n"
"FILES\n"
"    <Matrix> ................ I The matrix\n"
"    <Result> ................ O The transformed matrix\n"
};

static MtxApplication_t *App = NULL;
static const char *iname, *oname;
static MtxFile_t *ifile, *ofile;
static PTR m1;				/* Workspace, one row */



static int Init(int argc, const char **argv)

{
    /* Parse command line
       ------------------ */
    if ((App = AppAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    if (AppGetArguments(App,2,2) < 0)
	return -1;
    iname = App->ArgV[0];
    oname = App->ArgV[1];
    return 0;
}




static int OpenFiles()

{
    /* Open the input file
       ------------------- */
    if ((ifile = MfOpen(iname)) == NULL)
	return -1;
    if (ifile->Field < 2) 
    {
	MTX_ERROR2("%s: %E",iname,MTX_ERR_NOTMATRIX);
	return -1;
    }
    FfSetField(ifile->Field); 
    FfSetNoc(ifile->Noc);
    MESSAGE(0,("Characteristic is %d\n",FfChar));

    /* Open output file, allocate memory.
       ---------------------------------- */
    if ((ofile = MfCreate(oname,ifile->Field,ifile->Nor,ifile->Noc)) == NULL)
    {
	MTX_ERROR("Cannot create output file");
	return -1;
    }
    m1 = FfAlloc(1);
    if (m1 == NULL)
	return -1;
    return 0;
}



static void Cleanup()

{
    if (ifile != NULL)
	MfClose(ifile);
    if (ofile != NULL)
	MfClose(ofile);
    SysFree(m1);
    if (App != NULL)
	AppFree(App);
}


static int FrobeniusMap()

{
    int i;

    /* Apply the frobenius map to each entry.
       -------------------------------------- */
    for (i = 0; i < ifile->Nor; ++i)
    {
	int k;
	if (MfReadRows(ifile,m1,1) != 1)
	{
	    MTX_ERROR1("Error reading %s",iname);
	    return 1;
	}
	for (k = 0; k < FfNoc; ++k)
	{
	    FEL f1 = FfExtract(m1,k);
	    FEL f2 = f1;
	    int n;
	    for (n = FfChar - 1; n != 0; --n)
		f2 = FfMul(f1,f2);
	    FfInsert(m1,k,f2);
	}
	if (MfWriteRows(ofile,m1,1) != 1)
	{
	    MTX_ERROR1("Error writing %s",oname);
	    return 1;
	}
    }

    return 0;
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
    if (OpenFiles() != 0)
	return 1;
    if (FrobeniusMap() != 0)
	return 1;
    Cleanup();
    return 0;
}



/**
@page prog_zfr zfr - Frobenius Automorphism

@section zfr_syntax Command Line
<pre>
zfr @em Options @em Mat @em Result
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par @em Mat
  Input matrix
@par @em Result
  Result matrix

@section zfr_inp Input Files
@par @em Mat
  Input matrix

@section zfr_out Output Files
@par @em Result
  Result matrix

@section zfr_desc Description
This program reads a matrix, applies the Frobenius automorphism 
x↦x<sup>p</sup>, where p is the characteristic of the field, 
to each entry and writes out the result. 
*/
