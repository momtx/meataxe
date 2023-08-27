#include "meataxe.h"



static MtxApplicationInfo_t AppInfo = { 
"zcl", "Clean Matrix", 
"SYNTAX\n"
"    zcl <Subsp> <Mat> <Cleaned mat> <Ops>\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"\n"
"FILES\n"
"    <Subsp> ................. I The subspace to clean with\n"
"    <Mat> ................... I The matrix to be cleaned\n"
"    <Cleaned> ............... O The cleaned matrix\n"
"    <Ops> ................... O Row operations that were performed\n"
};
static MtxApplication_t *App = NULL;


const char *SpcName = NULL, *MatName = NULL, *OpName = NULL, *ClName = NULL;
static Matrix_t *Space;
static MtxFile_t *Matrix, *Cleaned, *Op;


static int Clean1()

{
    PTR op, row;
    int j;
    int rc = 0;

    /* Check compatibility.
       -------------------- */
    if (Matrix->Noc != Space->Noc || Matrix->Field != Space->Field)
    {
	mtxAbort(MTX_HERE,"%s and %s: %s",SpcName,MatName,MTX_ERR_INCOMPAT);
	return 1;
    }

    /* Build the pivot table.
       ---------------------- */
    if (matPivotize(Space) != 0)
    {
	mtxAbort(MTX_HERE,"%s: %s",SpcName,MTX_ERR_NOTECH);
	return -1;
    }

    /* Allocate work space.
       -------------------- */
    if ((op = ffAlloc(1,Space->Nor)) == NULL)
	rc = -1;
    if ((row = ffAlloc(1, Space->Noc)) == NULL)
	rc = -1;

    /* Process the matrix row by row.
       ------------------------------ */
    for (j = 0; rc == 0 && j < Matrix->Nor; ++j)
    {
	if (mfReadRows(Matrix,row,1) != 1)
	    rc = -1;
	ffCleanRow2(row,Space->Data,Space->Nor,Space->Noc,Space->PivotTable,op);
	if (mfWriteRows(Cleaned,row,1) != 1)
	    rc = -1;
	if (mfWriteRows(Op,op,1) != 1)
	    rc = -1;
	ffMulRow(op,FF_ZERO, Space->Noc);
    }

    /* Clean up.
       --------- */
    if (op != NULL)
	sysFree(op);
    if (row != NULL)
	sysFree(row);

    return rc;
}



int main(int argc, char **argv)

{
    int rc;


    if ((App = appAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    if (appGetArguments(App,4,4) < 0)
	return -1;
    SpcName = App->ArgV[0];
    MatName = App->ArgV[1];
    ClName = App->ArgV[2];
    OpName = App->ArgV[3];

    if ((Space = matLoad(SpcName)) == NULL)
	rc = -1;
    if ((Matrix = mfOpen(MatName)) == NULL)
	rc = -1;
    if ((Cleaned = mfCreate(ClName,ffOrder,Matrix->Nor,Matrix->Noc)) == NULL)
	rc = -1;
    if ((Op = mfCreate(OpName,ffOrder,Matrix->Nor,Space->Nor)) == NULL)
	rc = -1;

    /* Clean.
       ------ */
    rc = Clean1();

    /* Clean up.
       --------- */
    if (Space != NULL) matFree(Space);
    if (Matrix != NULL) mfClose(Matrix);
    if (Cleaned != NULL) mfClose(Cleaned);
    if (Op != NULL) mfClose(Op);
    if (App != NULL) appFree(App);

    return rc;
}


/**
@page prog_zcl zcl - Clean

@section zcl_syntax Command Line
<pre>
zcl [@em Options] @em Subsp @em Mat @em CleanedMat @em Ops
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par @em Subsp
  Subspace to clean with.
@par @em Mat
  Matrix to be cleaned.
@par @em CleanedMat
  Cleaned Matrix.
@par @em Ops
  Row operations.

@section ifile Input Files
@par @em Subsp
  Subspace to clean with, a N⨯M matrix in echelon form.
@par @em Mat
  Matrix to be cleaned (L⨯M).

@section ofile Output Files
@par @em CleanedMat
  Cleaned matrix (L⨯M).
@par @em Ops
  Row operations (L⨯N).


@section zcl_desc Description
This program "cleans" @em Mat with @em Subsp, i.e., it adds suitable linear
combinations of rows of @em Subsp to each row of @em Mat such that all pivot
columns in the result are zero. It writes two matrices to @em CleanedMat and @em Ops
such that @em Mat = @em Ops⋅@em Subsp + @em CleanedMat.
If @em Subsp is not in echelon form, it is
first reduced to echelon form, and the previous equation holds for the
reduced Matrix.
@em Subsp and @em Mat must be over the same field and have the same number of columns.

One use of this program is to calculate the action of a generator on an invariant
subspace: Take the subspace in echelon form (as it is on output from @ref prog_zsp "zsp"),
and multiply it by a generator. Cleaning the result with the original basis yields a
zero matrix, and @em RowOps is the action of the generator on the invariant subspace.
For example, if @c subsp is the subspace and @c gen is the generator,
<pre>
# zmu subsp gen image
# zcl subsp image null gen_s
</pre>
calculates the action on the subspace in @c gen_s.

The action on the quotient of a given subspace can be calculated with
the @ref prog_zqt "zqt" program.

@section imspl Implementation Details
The subspace is loaded into memory and, if necessary, reduced to echelon form.
The second matrix, is then processed a row at a time.
Row operations are performed to clear out the pivot points of the input row using the 
rows of the first input matrix. A row describing what was done is written 
out to @em Ops, and the remnant (clean) row is output to @em CleanedMat.
**/

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
