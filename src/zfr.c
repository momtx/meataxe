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
static MtxFile_t* ifile;
static MtxFile_t* ofile;


////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char **argv)
{
    App = appAlloc(&AppInfo,argc,argv);
    appGetArguments(App,2,2);
    iname = App->argV[0];
    oname = App->argV[1];
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void openFiles()
{
    ifile = mfOpen(iname);
    mfReadHeader(ifile);
    if (mfObjectType(ifile) != MTX_TYPE_MATRIX)
	mtxAbort(MTX_HERE,"%s: %s",iname,MTX_ERR_NOTMATRIX);
    ffSetField(ifile->header[0]); 
    MESSAGE(0,("Characteristic is %lu\n",(unsigned long) ffChar));
    ofile = mfCreate(oname,ifile->header[0],ifile->header[1],ifile->header[2]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void cleanup()
{
   mfClose(ifile);
   mfClose(ofile);
   appFree(App);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void frobeniusMap()
{
    const uint32_t nor = ifile->header[1];
    const uint32_t noc = ifile->header[2];
    PTR m1 = ffAlloc(1, noc);

    for (uint32_t i = 0; i < nor; ++i)
    {
	mfReadRows(ifile,m1,1, noc);
	for (uint32_t k = 0; k < noc; ++k)
	{
	    FEL f1 = ffExtract(m1,k);
	    FEL f2 = f1;
	    for (uint32_t n = ffChar - 1; n != 0; --n)
		f2 = ffMul(f1,f2);
	    ffInsert(m1,k,f2);
	}
	mfWriteRows(ofile,m1,1,noc);
    }
    sysFree(m1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    init(argc,argv);
    openFiles();
    frobeniusMap();
    cleanup();
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// *INDENT-OFF*

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
