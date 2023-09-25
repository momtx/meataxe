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

const char* SpcName = NULL;
const char* MatName = NULL;
const char* OpName = NULL;
const char* ClName = NULL;
static Matrix_t* Space = NULL;
static MtxFile_t* matrixFile = NULL;
static MtxFile_t* cleanedFile = NULL;
static MtxFile_t* opFile = NULL;
static uint32_t noc = 0;
static uint32_t spaceNor = 0;
static uint32_t matrixNor = 0;

////////////////////////////////////////////////////////////////////////////////////////////////////

static void clean()
{
   matPivotize(Space);
   PTR op = ffAlloc(1,spaceNor);
   PTR row = ffAlloc(1, noc);

   for (uint32_t j = 0; j < matrixNor; ++j)
   {
      mfReadRows(matrixFile, row, 1, noc);
      ffMulRow(op, FF_ZERO, spaceNor);
      ffCleanRow2(row,Space->data,spaceNor,noc,Space->pivotTable,op);
      mfWriteRows(cleanedFile,row,1, noc);
      mfWriteRows(opFile,op,1,spaceNor);
   }
   sysFree(op);
   sysFree(row);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char **argv)
{
    App = appAlloc(&AppInfo,argc,argv);
    appGetArguments(App,4,4);
    SpcName = App->argV[0];
    MatName = App->argV[1];
    ClName = App->argV[2];
    OpName = App->argV[3];

    Space = matLoad(SpcName);
    noc = Space->noc;
    spaceNor = Space->nor;

    matrixFile = mfOpen(MatName);
    mfReadHeader(matrixFile);
    if (matrixFile->header[0] != Space->field || matrixFile->header[2] != noc)
       mtxAbort(MTX_HERE, "%s and %s: %s", SpcName, MatName, MTX_ERR_INCOMPAT);
    matrixNor = matrixFile->header[1];

    cleanedFile = mfCreate(ClName,ffOrder,matrixNor,noc);
    opFile = mfCreate(OpName,ffOrder,matrixNor,spaceNor);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   init(argc, argv);
   clean();

   matFree(Space);
   mfClose(matrixFile);
   mfClose(cleanedFile);
   mfClose(opFile);
   appFree(App);

   return 0;
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
