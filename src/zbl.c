////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - 
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static MtxApplicationInfo_t AppInfo = {
    "zbl", "Bottom Left of a Matrix",
    "SYNTAX\n"
    "    zbl " MTX_COMMON_OPTIONS_SYNTAX " <Matrix> <Result>\n"
    "\n"
    "ARGUMENTS\n"
    "    <Matrix> ................ Input matrix\n"
    "    <Result> ................ Output matrix\n"
    "\n"
    "OPTIONS\n"
    MTX_COMMON_OPTIONS_DESCRIPTION
};

static MtxApplication_t *App = NULL;
static const char *iname = 0, *oname = 0;



static int Init(int argc, char **argv)
{
    if ((App = appAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;

    /* Parse command line.
       ------------------- */
    if (appGetArguments(App,2,2) < 0)
	return -1;
    
    iname = App->argV[0];
    oname = App->argV[1];
    return 0;
}



/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char *argv[])
{
    if (Init(argc, (char **)argv) != 0)
    {
	mtxAbort(MTX_HERE,"Initialization failed");
	return 1;
    }

    // Open files.
    MtxFile_t* inputFile = mfOpen(iname, "rb");
    mfReadHeader(inputFile);
    if (mfObjectType(inputFile) != MTX_TYPE_MATRIX)
	mtxAbort(MTX_HERE,"%s: %s", iname, MTX_ERR_NOTMATRIX);
    ffSetField(inputFile->header[0]);
    const uint32_t fieldOrder = inputFile->header[0];
    const uint32_t nor = inputFile->header[1];
    const uint32_t noc = inputFile->header[2];
    MtxFile_t* outputFile = mfCreate(oname, fieldOrder, nor, noc);

    // Do the job.
    PTR m1 = ffAlloc(1, noc);
    for (uint32_t i = 0; i < nor; ++i)
    {
	ffReadRows(inputFile, m1, 1, noc);
	for (uint32_t j = i + 1; j < noc; ++j)
	    ffInsert(m1,j,FF_ZERO);
	ffWriteRows(outputFile, m1, 1, noc);
    }
    sysFree(m1);

    mfClose(inputFile);
    mfClose(outputFile);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// *INDENT-OFF*

/**
@page prog_zbl zbl - Bottom Left

@section zbl_syntax Command Line
<pre>
zbl [@em Options] @em Input @em Output
</pre>

@par @em Options
Standard options, see @ref prog_stdopts

@par @em Input
Input matrix

@par @em Output
Result matrix

@section zbl_inp Input Files
@par @em Input
Input matrix

@section zbl_out Output Files
@par @em Output
Result matrix

@section zbl_desc Description

This program reads in a matrix, zeroizes all entries above the main diagonal,
and writes out the result. For example:
<pre>
        Mat                  Result
12100201201212212122  10000000000000000000
21000000000000001111  21000000000000000000
12212120101201201201  12200000000000000000
01020010020102012012  01020000000000000000
12100012101201012012  12100000000000000000
12120212121012012012  12120200000000000000
</pre>
The input matrix need not be square, and the output matrix has always the same dimensions as the
input matrix.

If file names are omitted, the matrix is read from G1 and output goes to P2.

The purpose of this program is to enable the MeatAxe to check if an irreducible representation
in characteristic 2 fixes a quadratic form. This job is not particularly simple -- in many ways
it is just a bodge, but it is possible:
-# Put the representation into standard basis using @ref prog_zsp "zsp" in "standard basis" mode.
-# Find the symplectic form fixed by using @ref prog_zsp "zsp" to make the matrix that conjugates
  the representation to its dual.
-# Quadratic forms can be represented by lower triangular matrices. Since the representation
  is in standard basis (so all the basis vectors are images of the first under the group) the
  diagonal entries of any fixed quadratic form must all be equal, so try each field entry in
  turn by adding that scalar matrix to the bottom half of the symplectic form.
-# For each quadratic form Q made as in 3, test if it is fixed by forming G<sup>T</sup>QG for each
  generator G, and checking that the diagonal is still the same as it was before (the
  symplectic form should have been checked before starting). The check can be done by adding
  the form to the result, then doing @b zbl, @b ztr and @b zbl again --- the result will be the zero
  matrix (use @ref prog_znu "znu") iff the form was fixed (given that the symplectic one was).

*/
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
