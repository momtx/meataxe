////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Apply the Frobenius automorphism to a matrix.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */


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



static int Init(int argc, char **argv)

{
    /* Parse command line
       ------------------ */
    if ((App = appAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    if (appGetArguments(App,2,2) < 0)
	return -1;
    iname = App->ArgV[0];
    oname = App->ArgV[1];
    return 0;
}




static int OpenFiles()
{
    /* Open the input file
       ------------------- */
    if ((ifile = mfOpen(iname)) == NULL)
	return -1;
    if (ifile->Field < 2) 
    {
	mtxAbort(MTX_HERE,"%s: %s",iname,MTX_ERR_NOTMATRIX);
	return -1;
    }
    ffSetField(ifile->Field); 
    ffSetNoc(ifile->Noc);
    MESSAGE(0,("Characteristic is %d\n",ffChar));

    /* Open output file, allocate memory.
       ---------------------------------- */
    if ((ofile = mfCreate(oname,ifile->Field,ifile->Nor,ifile->Noc)) == NULL)
    {
	mtxAbort(MTX_HERE,"Cannot create output file");
	return -1;
    }
    m1 = ffAlloc(1, ifile->Noc);
    if (m1 == NULL)
	return -1;
    return 0;
}



static void Cleanup()

{
    if (ifile != NULL)
	mfClose(ifile);
    if (ofile != NULL)
	mfClose(ofile);
    sysFree(m1);
    if (App != NULL)
	appFree(App);
}


static int FrobeniusMap()

{
    int i;

    /* Apply the frobenius map to each entry.
       -------------------------------------- */
    for (i = 0; i < ifile->Nor; ++i)
    {
	int k;
	if (mfReadRows(ifile,m1,1) != 1)
	{
	    mtxAbort(MTX_HERE,"Error reading %s",iname);
	    return 1;
	}
	for (k = 0; k < ffNoc; ++k)
	{
	    FEL f1 = ffExtract(m1,k);
	    FEL f2 = f1;
	    int n;
	    for (n = ffChar - 1; n != 0; --n)
		f2 = ffMul(f1,f2);
	    ffInsert(m1,k,f2);
	}
	if (mfWriteRows(ofile,m1,1) != 1)
	{
	    mtxAbort(MTX_HERE,"Error writing %s",oname);
	    return 1;
	}
    }

    return 0;
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
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
