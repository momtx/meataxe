#include <meataxe.h>


MTX_DEFINE_FILE_INFO

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
	MTX_ERROR3("%s and %s: %E",SpcName,MatName,MTX_ERR_INCOMPAT);
	return 1;
    }

    /* Build the pivot table.
       ---------------------- */
    if (MatPivotize(Space) != 0)
    {
	MTX_ERROR2("%s: %E",SpcName,MTX_ERR_NOTECH);
	return -1;
    }

    /* Allocate work space.
       -------------------- */
    FfSetNoc(Space->Nor);
    if ((op = FfAlloc(1)) == NULL)
	rc = -1;
    FfSetNoc(Space->Noc);
    if ((row = FfAlloc(1)) == NULL)
	rc = -1;

    /* Process the matrix row by row.
       ------------------------------ */
    for (j = 0; rc == 0 && j < Matrix->Nor; ++j)
    {
	if (MfReadRows(Matrix,row,1) != 1)
	    rc = -1;
	FfCleanRow2(row,Space->Data,Space->Nor,Space->PivotTable,op);
	if (MfWriteRows(Cleaned,row,1) != 1)
	    rc = -1;
	if (MfWriteRows(Op,op,1) != 1)
	    rc = -1;
	FfMulRow(op,FF_ZERO);
    }

    /* Clean up.
       --------- */
    if (op != NULL)
	SysFree(op);
    if (row != NULL)
	SysFree(row);

    return rc;
}



int main(int argc, const char **argv)

{
    int rc;


    if ((App = AppAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    if (AppGetArguments(App,4,4) < 0)
	return -1;
    SpcName = App->ArgV[0];
    MatName = App->ArgV[1];
    ClName = App->ArgV[2];
    OpName = App->ArgV[3];

    if ((Space = MatLoad(SpcName)) == NULL)
	rc = -1;
    if ((Matrix = MfOpen(MatName)) == NULL)
	rc = -1;
    if ((Cleaned = MfCreate(ClName,FfOrder,Matrix->Nor,Matrix->Noc)) == NULL)
	rc = -1;
    if ((Op = MfCreate(OpName,FfOrder,Matrix->Nor,Space->Nor)) == NULL)
	rc = -1;

    /* Clean.
       ------ */
    rc = Clean1();

    /* Clean up.
       --------- */
    if (Space != NULL) MatFree(Space);
    if (Matrix != NULL) MfClose(Matrix);
    if (Cleaned != NULL) MfClose(Cleaned);
    if (Op != NULL) MfClose(Op);
    if (App != NULL) AppFree(App);

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

