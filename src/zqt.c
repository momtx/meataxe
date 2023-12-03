////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Projection on quotient.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>
#include <stdlib.h>




/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */


static MtxApplicationInfo_t AppInfo = { 
"zqt", "Clean And Quotient", 
"SYNTAX\n"
"    zqt [-i] [<Subsp> <Matrix> <Quotient>]\n"
"\n"
"OPTIONS\n"
"    -i   Take only insignificant rows of <Matrix>. <Quotient> will be\n"
"         the action of <Matrix> on the quotient by <Subspace>.\n"
"\n"
"FILES\n"
"    <Subsp>    i  The invariant subspace, in semi-echelon form\n"
"    <Matrix>   i  The matrix, must have the same number of columns\n"
"    <Quotient> o  Insignificant columns of <Matrix>, cleaned with <Subspace>\n"
"\n"
};


static MtxApplication_t *App = NULL;
static int opt_i = 0;

static const char *fileNameS;
static Matrix_t *S = NULL;

static const char *fileNameM;
static MtxFile_t *fileM = NULL;
static int32_t norM = 0;
static int32_t nocM = 0;
static PTR bufferM = NULL;

static const char *fileNameQ;
static MtxFile_t *fileQ = NULL;
static PTR bufferQ = NULL;

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char **argv)
{
    App = appAlloc(&AppInfo,argc,argv);
    opt_i = appGetOption(App,"-i");
    appGetArguments(App,3,3);
    fileNameS = App->argV[0];
    fileNameM = App->argV[1];
    fileNameQ = App->argV[2];
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int readFiles()
{
   // subspace
    S = matLoad(fileNameS);
    if (S == NULL)
	return 1;
    matPivotize(S);

    // matrix
    fileM = mfOpen(fileNameM, "rb");
    mfReadHeader(fileM);
    if (mfObjectType(fileM) != MTX_TYPE_MATRIX)
       mtxAbort(MTX_HERE, "%s: %s", fileNameM, MTX_ERR_NOTMATRIX);
    norM = fileM->header[1];
    nocM = fileM->header[2];
    if (nocM != S->noc)
	mtxAbort(MTX_HERE,"%s and %s: %s",fileNameS,fileNameM,MTX_ERR_INCOMPAT);
    if (opt_i && norM != nocM)
       mtxAbort(MTX_HERE, "%s: %s", fileNameM, MTX_ERR_NOTSQUARE);
    
    bufferM = ffAlloc(1, nocM);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int isPivot(int i)
{
   uint32_t *piv = S->pivotTable;
   uint32_t dimS = S->nor;
   for (uint32_t k = 0; k < dimS; ++k)
   {
      if (i == piv[k])
         return 1;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void doit()
{
    int i;
    const uint32_t *non_pivot = S->pivotTable + S->nor;
    const uint32_t quotientDim = S->noc - S->nor;

    const uint32_t norM = fileM->header[1];
    const uint32_t nocM = fileM->header[2];

    // output file
    fileQ = mfCreate(fileNameQ,ffOrder,opt_i ? quotientDim : norM, quotientDim);
    bufferQ = ffAlloc(1, nocM);

    for (i = 0; i < norM; ++i)
    {
	ffReadRows(fileM,bufferM, 1, nocM);

	// When calculating the action, take only insignificant rows.
	if (opt_i && isPivot(i))
	    continue;

	// Clean and extract insignificant columns.
	ffCleanRow(bufferM,S->data,S->nor,S->noc,S->pivotTable);
	ffMulRow(bufferQ,FF_ZERO, quotientDim);
	for (uint32_t k = 0; k < quotientDim; ++k)
	    ffInsert(bufferQ,k,ffExtract(bufferM,non_pivot[k]));

	// Write the output row.
	ffWriteRows(fileQ,bufferQ, 1, quotientDim);
    }

    sysFree(bufferQ);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void cleanup()
{
   matFree(S);
   mfClose(fileM);
   mfClose(fileQ);
   appFree(App);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    init(argc,argv);
    readFiles();
    doit();
    cleanup();
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// *INDENT-OFF*

/**
@page prog_zqt zqt - Clean and Quotient

@section zqt_syntax Command Line
<pre>
zqt [@em Options] [-i] @em Subsp @em Matrix @em Quot
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par -i
  Calculate the action on the quotient.
@par @em Subsp
  The subspace.
@par @em Mat
  The matrix.
@par @em Quot.
  The quotient.

@section zqt_inp Input Files
@par @em Subsp
  The subspace, a M×N matrix in echelon form.
@par @em Mat
  The matrix (L×N).

@section zqt_out Output Files
@par @em Quot.
  The quotient, a L×(N-M) matrix.
  With "-i", the action on the quotient, a (N-M)×(N-M) matrix.

@section zqt_desc Description
This program reads in a subspace and applies the canonical map to its quotient to a matrix.
The result is written out to @em Quot.
@em Subsp must be a matrix in semi-echelon form and have the same number of columns as
@em Mat.

Otherwise the program reads in @em Subsp, builds a table of pivot
columns and then proceeds, row by row, through @em Matrix. For
each row, the significant entries are zeroized by adding the correct
multiple of rows of @em Subsp. The insignificant columns are then
extracted and written out to @em Quot. Hence
- @em Subsp has M rows, N columns and is in echelon form,
- @em Matrix has L rows, N columns and is otherwise arbitrary, and
- @em Quot has L rows and N-M columns.
In other words, the program calculates the projection of @em Matrix
onto the B-A dimensional quotient space defined by @em Subsp.  
If the "-i" option is used, @b zqt calculates the action of the
Matrix on the quotient. This is done by projecting the matrix as
explained above, and taking only the @em insignificant rows.
Insignificant rows are defined by treating
the pivot table as a table of rows rather than columns.
Example: Let "spc" be an invariant subspace and "z1" an
algebra element (a square matrix). Then, after
<pre>
zqt -i spc z1 q1</pre>
"q1" contains the action of "z1" on the quotient by "spc".

Another, less obvious use of @b zqt is to condense a matrix representation.
First, find an element E of the group algebra with stable rank, i.e.,
rank(E*E) = rank(E). This can be done by taking any element F of the group
algebra and raising it to higher powers until the rank stabilizes.
We may then condense onto the kernel of E as follows
<pre>
zef E X        X is the echelon form of Image(E)
znu E Y        Y is the kernel of E
zqt X Y Z      calculate the canonical projection of Y ...
ziv Z T        ... and adjust Y so that the canonical ...
zmu T Y Y1     ... projection of Y1 is the identity
zmu Y1 Z1 T1   calculate KZ1 = condensed Z1
zqt X T1 KZ1
zmu Y1 Z2 T1   calculate KZ2 = condensed Z2
zqt X T1 KZ2
</pre>


@section zqt_impl Implementation Details
It is not completely checked that @em Subsp is in echelon form.

The Subspace and one row of both @em Matrix and @em Subsp must fit into memory.
**/

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
