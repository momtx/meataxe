////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Multiply matrices or permutations.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>




/* --------------------------------------------------------------------------
   Global data
   -------------------------------------------------------------------------- */


static MtxApplicationInfo_t AppInfo = { 
"zmu", "Multiply", 
"SYNTAX\n"
"    zmu <A> <B> <Result>\n"
"\n"
"FILES\n"
"    <A> and <B> are the objects to be multiplied. Their product\n"
"    (A*B) is written to <Result>. Compatible data types are:\n"
"\n"
"        M(a,b) * M(b,c)                   = M(a,c)\n"
"        M(1,1) * M(a,b) = M(a,b) * M(1,1) = M(a,b)\n"
"        P(a) * P(b)                       = P(max {a,b})\n"
"        M(a,b) * P(b)                     = M(a,b)\n"
"        P(a) * M(a,b)                     = M(a,b)\n"
"\n"
"    M(a,b) means `a by b matrix' and P(a) `Permutation of degree a'\n"
};


static MtxApplication_t* App = NULL;
static const char* fileNameA;
static const char* fileNameB;
static const char* fileNameC;
static MtxFile_t* fileA = NULL; 
static MtxFile_t* fileB = NULL;
static MtxFile_t* fileC = NULL;

////////////////////////////////////////////////////////////////////////////////////////////////////

// Multiply permutation * matrix

static int multpm(void)
{
    const uint32_t degreeA = fileA->header[1];
    const uint32_t fieldB = fileB->header[0];
    const uint32_t norB = fileB->header[1];
    const uint32_t nocB = fileB->header[2];
    if (degreeA != norB) 
	mtxAbort(MTX_HERE,"%s and %s: %s",fileNameA,fileNameB,MTX_ERR_INCOMPAT);

    // read input files
    Perm_t* permA = permReadData(fileA->file, fileA->header);
    Matrix_t *matrixB = matAlloc(fieldB, norB, nocB);

    // Write out the rows of <B> in the order defined by <A>.
    fileC = mfCreate(fileNameC,fieldB, norB, nocB);
    for (uint32_t i = 0; i < degreeA; ++i)
    {	
       PTR row = ffGetPtr(matrixB->data, permA->data[i], nocB);
       mfWriteRows(fileC, row, 1, nocB);
    }

    permFree(permA);
    matFree(matrixB);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Multiply matrix * permutation

static void multmp(void)
{
   const uint32_t fieldA = fileA->header[0];
   const uint32_t norA = fileA->header[1];
   const uint32_t nocA = fileA->header[2];
   const uint32_t degreeB = fileB->header[1];
   if (nocA != degreeB) 
      mtxAbort(MTX_HERE,"%s and %s: %s", fileNameA, fileNameB, MTX_ERR_INCOMPAT);

   // Read the permutation (B)
   Perm_t* perm = permReadData(fileB->file, fileB->header);

   // Allocate workspace (two rows of A).
   ffSetField(fieldA);
   PTR row_in = ffAlloc(2, nocA);
   PTR row_out = ffGetPtr(row_in, 1, nocA);

   // Create the output file.
   fileC = mfCreate(fileNameC,fieldA, norA, nocA);

   // Process A row by row. Permute the marks of each row according to B.
   for (uint32_t i = 0; i < norA; ++i)
   {
      mfReadRows(fileA, row_in, 1, nocA);
      ffPermRow(row_in, perm->data, nocA, row_out);
      mfWriteRows(fileC, row_out,1, nocA);
   }

   permFree(perm);
   ffFree(row_in);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

// Multiply scalar with matrix.

static int multsm(MtxFile_t* fileS, MtxFile_t* fileM)
{
    const uint32_t field = fileS->header[0];
    ffSetField(field);
    PTR ms = ffAlloc(1,1);
    mfReadRows(fileS,ms,1,1);
    const FEL f = ffExtract(ms,0);
    ffFree(ms);

    const uint32_t norM = fileM->header[1];
    const uint32_t nocM = fileM->header[2];
    PTR mm = ffAlloc(1, nocM);
	
    fileC = mfCreate(fileNameC,field ,norM,nocM);
    for (int i = 0; i < norM; ++i)
    {	
	mfReadRows(fileM, mm, 1, nocM);
	ffMulRow(mm, f, nocM);
	mfWriteRows(fileC, mm,1, nocM);
    }
    sysFree(mm);
    return 0;
}



/* ------------------------------------------------------------------
   multmm() - Multiply two matrices
   ------------------------------------------------------------------ */

static void multmm(void)
{
    uint32_t fieldA = fileA->header[0];
    if (fileB->header[0] != fieldA)
	mtxAbort(MTX_HERE,"%s and %s: %s (different fields)",fileNameA,fileNameB, MTX_ERR_INCOMPAT);

    const uint32_t norA =  fileA->header[1];
    const uint32_t nocA =  fileA->header[2];
    const uint32_t norB =  fileB->header[1];
    const uint32_t nocB =  fileB->header[2];

    if (norA == 1 && nocA == 1) {
	multsm(fileA,fileB);
        return;
    }
    if (norB == 1 && nocB == 1) {
	multsm(fileB,fileA);
        return;
    }
    if (nocA != norB) 
	mtxAbort(MTX_HERE,"%s and %s: %s",fileNameA,fileNameB,MTX_ERR_INCOMPAT);
   
    ffSetField(fieldA);
    PTR rowA = ffAlloc(1, nocA);
    PTR matrixB = ffAlloc(norB, nocB);
    mfReadRows(fileB,matrixB, norB, nocB);
    PTR rowC = ffAlloc(1, nocB);

    fileC = mfCreate(fileNameC, fieldA, norA, nocB);
    for (uint32_t i = 0; i < norA; ++i)
    {
        mfReadRows(fileA, rowA, 1, nocA);
	ffMapRow(rowA, matrixB, norB, nocB, rowC);
	mfWriteRows(fileC,rowC, 1, nocB);
    }
    sysFree(rowC);
    sysFree(matrixB);
    sysFree(rowA);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Multiply two permutations

static void multpp(void)
{
    Perm_t *permA = permReadData(fileA->file, fileA->header);
    Perm_t *permB = permReadData(fileB->file, fileB->header);
          
    permMul(permA,permB);
    permSave(permA,fileNameC);

    permFree(permA);
    permFree(permB);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char **argv)
{
    App = appAlloc(&AppInfo,argc,argv);
    appGetArguments(App,3,3);

    fileNameA = App->argV[0];
    fileNameB = App->argV[1];
    fileNameC = App->argV[2];
    if (!strcmp(fileNameA,fileNameC) || !strcmp(fileNameB,fileNameC))
	mtxAbort(MTX_HERE,"Output file would overwrite input file");

    // Open input files
    fileA = mfOpen(fileNameA);
    mfReadHeader(fileA);
    fileB = mfOpen(fileNameB);
    mfReadHeader(fileB);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void Cleanup()
{
   mfClose(fileA);
   mfClose(fileB);
   if (fileC != NULL) 
      mfClose(fileC);
   appFree(App);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{   
    init(argc,argv);

    // Call the appropriate multiplication function.
    const uint32_t typeA = mfObjectType(fileA);
    const uint32_t typeB = mfObjectType(fileB);
    if (typeA == MTX_TYPE_PERMUTATION && typeB == typeA)
	multpp(); 
    else if (typeA == MTX_TYPE_MATRIX && typeB == typeA)
	multmm(); 
    else if (typeA == MTX_TYPE_MATRIX && typeB == MTX_TYPE_PERMUTATION) 
	multmp(); 
    else if (typeA == MTX_TYPE_PERMUTATION && typeB == MTX_TYPE_MATRIX)
	multpm(); 
    else
	mtxAbort(MTX_HERE,"%s and %s: %s",fileNameA,fileNameB,MTX_ERR_INCOMPAT);

    Cleanup();

    return 0;
}




/**
@page prog_zmu zmu - Multiply

@see  @ref prog_zpt

@section zmu_syntax Command Line
<pre>
zmu @em [Options] [-r @em Row[/@em NRows]] [-c @em Col[/@em NCols]] @em A @em B @em Result
</pre>

@par @em Options
Standard options, see @ref prog_stdopts

@par -r @em Row[/@em NRows]
Divide the matrix @em A horizontally into @em NRows slices (default: 2) and use the
@em Row-th slice as the left factor. @em NRows must not be larger than the number of rows
of @em A.

@par -c @em Col[/@em NCols]
Divide the matrix @em B vertically into @em NCols slices (default: 2) and use the
@em Col-th slice as right factor. @em NCols must not be larger than the number of 
columns of @em B.

@par @em A
Left factor.

@par @em B
Right factor.

@par @em Result
Product.

@section zmu_inp Input Files

@par @em A
Left factor.

@par @em B
Right factor.

@section zmu_out Output Files

@par @em Result
Product.

@section zmu_desc Description
This program reads two matrices or permutations and writes their product to @em Result.

The input files must contain two compatible objects, i.e., their product must be defined.
Currently, @b zmu can handle the following data types:

- Both files are matrices over the same field, and the number of columns of @em A
  equals the number of rows of @em B. In this case, @b zmu calculates the standard
  matrix product.

- One of the operands is a one by one matrix, and the other is
  any matrix over the same field. In this case, the one by one
  matrix is interpreted as a scalar, and the program calculates
  the corresponding multiple of the matrix.

- Both input files are permutations of degree a and b, respectively.
  The result is a permutation C of degree max(a,b), which is defined
  by C(x) = B(A(x)). If the permutations are of different degrees,
  the smaller permutation is extended to the larger degree by adding fixed points.

- @em A is a matrix, @em B is a permutation and the degree of the permutation equals the
  number of columns of the matrix.
  The result is a matrix of the same size which is calculated
  from the input matrix by permuting the marks of each row in
  the following way: The i-th mark of the row is stored as
  the k-th mark of the result if the permutation maps i to k.

- @em A is a permutation of degree m, and @em B is a m by n matrix.
  The result is again a m by n matrix which consists of the rows of the input matrix,
  rearranged according to the permutation. If the permutation maps i to k, then the k-th
  row of the input matrix becomes the i-th row of the output matrix.
  Here is an example:
<pre>
            | 1 1 |     | 2 2 |
(1 2 3)  *  | 2 2 |  =  | 3 3 |
            | 3 3 |     | 1 1 |
</pre>

With these conventions, products between matrices and permutations are defined in a consistent way.
The associative law a(bc)=(ab)c holds whenever ab and bc are defined (a,b,c being matrices or
permutations). A permutation matrix created with @ref prog_zcv "zcv" or @ref prog_zcf "zcf",
if multiplied with another matrix, produces the same result as the original permutation.

@subsection bl Blockwise Matrix Multiplication
In the case of two matrices, a blockwise multiplication can be performed using the
"-r" and "-c" options. If one or both of these options are specified on the command line,
@em zmu will read only some rows of @em A and/or some columns of @em B.
Multiplying the two pieces together yields a rectangular piece of the
result. By default the result is divided into 4 pieces of (almost)
equal size. To calculate the 4 pieces successively, type
<pre>
zmu -r 1 -c 1 m1 m2 tmp11
zmu -r 1 -c 2 m1 m2 tmp12
zmu -r 2 -c 1 m1 m2 tmp21
zmu -r 2 -c 2 m1 m2 tmp22
</pre>
The resulting matrices `tmpXX' can then be pasted together using @ref prog_zpt "zpt":
<pre>
zpt -R 2 -C 2 result tmp
</pre>
This procedure can be used in a multi-processor environment
where each piece of the result is computed on a separate machine.

By adding an additional parameter to "-r" and/or "-c" you can
control the number of vertical or horizontal slices. For example,
<pre>
zmu -r 3/5
</pre>
means to cut @em A horizontally into five slices and use the third
slice for multiplication. The number of slice must not be greater
than the number of rows.
*/
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
