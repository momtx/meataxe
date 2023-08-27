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
static const char *mname, *sname, *oname;
static Matrix_t *Subspace = NULL;
static MtxFile_t *InputFile = NULL;
static MtxFile_t *OutputFile = NULL;
static PTR InputBuffer = NULL;
static PTR OutputBuffer = NULL;
static int QuotientDim;



/* ------------------------------------------------------------------
   init() - Process command line options and arguments
   ------------------------------------------------------------------ */


static int Init(int argc, char **argv)

{
    App = appAlloc(&AppInfo,argc,argv);
    if (App == NULL)
	return -1;
    opt_i = appGetOption(App,"-i");
    if (appGetArguments(App,3,3) < 0)
	return -1;
    sname = App->ArgV[0];
    mname = App->ArgV[1];
    oname = App->ArgV[2];
    return 0;
}




/* ------------------------------------------------------------------
   readfiles() - Read the input files, allocate tables
   ------------------------------------------------------------------ */

static int ReadFiles()

{
    /* Read the subspace, and build the pivot table.
       --------------------------------------------- */
    Subspace = matLoad(sname);
    if (Subspace == NULL)
	return 1;
    if (matPivotize(Subspace) < 0)
    {
	mtxAbort(MTX_HERE,"%s: %s",mname,MTX_ERR_NOTECH);
	return 1;
    }

    /* Open the input file, check compatibility, and allocate buffer.
       -------------------------------------------------------------- */
    InputFile = mfOpen(mname);
    if (InputFile == NULL)
	return 1;
    if (InputFile->Field != Subspace->Field || InputFile->Noc != Subspace->Noc)
    {
	mtxAbort(MTX_HERE,"%s and %S: %s",sname,mname,MTX_ERR_INCOMPAT);
	return 1;
    }
    if (opt_i && InputFile->Nor != InputFile->Noc) 
    {
	mtxAbort(MTX_HERE,"%s: %s",mname,MTX_ERR_NOTSQUARE);
	return 1;
    }
    InputBuffer = ffAlloc(1, InputFile->Noc);
    QuotientDim = Subspace->Noc - Subspace->Nor;

    /* Open the output file and allocate buffer.
       ----------------------------------------- */
    OutputFile = mfCreate(oname,ffOrder,opt_i ? QuotientDim : InputFile->Nor,
	QuotientDim);
    OutputBuffer = ffAlloc(1, InputFile->Noc);

    return 0;
}





/* ------------------------------------------------------------------
   IsPivot() - Find out if a given column is significant or not.
   ------------------------------------------------------------------ */

static int IsPivot(int i)

{
    int k;
    int *piv = Subspace->PivotTable;
    int sdim = Subspace->Nor;
    for (k = 0; k < sdim; ++k)
    {
	if (i == piv[k])
	    return 1;
    }
    return 0;
}




/* ------------------------------------------------------------------
   doit()
   ------------------------------------------------------------------ */

static int doit()

{
    int i;
    int *non_pivot = Subspace->PivotTable + Subspace->Nor;

    for (i = 0; i < InputFile->Nor; ++i)
    {
	int k;

	/* Read one row from the input file.
	   --------------------------------- */
	if (mfReadRows(InputFile,InputBuffer,1) != 1)
	{
	    mtxAbort(MTX_HERE,"Error reading vector");
	    return 1;
	}

	/* When calculating the action, take only `insignifican' rows.
	   ----------------------------------------------------------- */
	if (opt_i && IsPivot(i))
	    continue;

	/* Clean and extract insignificant columns.
	   ---------------------------------------- */
	ffCleanRow(InputBuffer,Subspace->Data,Subspace->Nor,Subspace->Noc,Subspace->PivotTable);
	ffMulRow(OutputBuffer,FF_ZERO, QuotientDim);
	for (k = 0; k < QuotientDim; ++k)
	    ffInsert(OutputBuffer,k,ffExtract(InputBuffer,non_pivot[k]));

	/* Write the output row.
	   --------------------- */
	if (mfWriteRows(OutputFile,OutputBuffer,1) != 1)
	{
	    mtxAbort(MTX_HERE,"Error writing vector");
	    return 1;
	}

    }

    return 0;
}



static void Cleanup()

{
    if (Subspace != NULL)
	matFree(Subspace);
    if (InputFile != NULL)
	mfClose(InputFile);
    if (OutputFile != NULL)
	mfClose(OutputFile);
    appFree(App);
}


/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char **argv)

{
    int rc;
    if (Init(argc,argv) != 0)
    {
	mtxAbort(MTX_HERE,"Initialization failed");
	return 1;
    }
    if (ReadFiles() != 0)
    {
	mtxAbort(MTX_HERE,"Error reading files");
	return 1;
    }
    rc = doit();
    Cleanup();
    return rc;
}


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
This program reads in a subspace and applies the canonical map to its quotient
on a matrix. The result is written out to @em Quot.
@em Subsp should be a matrix in semi-echelon form, and the two
input matrices must have the same field parameter and the same number
of columns. If this is not the case the program stops with an error message.

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
