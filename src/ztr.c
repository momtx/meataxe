////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Transpose a matrix.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"



/* ------------------------------------------------------------------
   Global Data
   ------------------------------------------------------------------ */


static const char *fileNameIn, *fileNameOut;


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

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    App = appAlloc(&AppInfo,argc,argv);
    appGetArguments(App,2,2);
    fileNameIn = App->argV[0];
    fileNameOut = App->argV[1];

    Matrix_t* matrixInp = matLoad(fileNameIn);
    const uint32_t field = matrixInp->field;
    const uint32_t norIn = matrixInp->nor;
    const uint32_t nocIn = matrixInp->noc;
    const PTR matrixIn = matrixInp->data;

    PTR rowOut = ffAlloc(1, norIn);
    MtxFile_t* fileOut = mfCreate(fileNameOut, field, nocIn, norIn);
    for (uint32_t j = 0; j < nocIn; ++j)
    {
	ffMulRow(rowOut, FF_ZERO, norIn);
	ffExtractColumn(matrixIn, norIn, nocIn, j, rowOut);
	mfWriteRows(fileOut, rowOut, 1, norIn);
    }
    mfClose(fileOut);
    sysFree(rowOut);
    matFree(matrixInp);
    appFree(App);
    return 0;
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
