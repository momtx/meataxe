////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Map vector under tensor product.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"
#include <stdlib.h>
#include <string.h>


/* ------------------------------------------------------------------
   Variables
   ------------------------------------------------------------------ */


static MtxApplicationInfo_t AppInfo = { 
"ztm", "Tensor Multiply",
"SYNTAX\n"
"    ztm [-QV] [-T <MaxTime>] <Vectors> <A> <B> <Result>"
"\n"
"ARGUMENTS\n"
"    <Vectors> ............... Input file: Vectors\n"
"    <A> ..................... Input file: Left factor (square matrix)\n"
"    <B> ..................... Input file: Right factor (square matrix)\n"
"    <Result> ................ Output file: Result\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
};
static MtxApplication_t *App = NULL;
static const char *fileNameA, *fileNameB;	/* Matrices */
static const char *fileNameVin, *fileNameVout;	/* Vectors, Result */
static Matrix_t *matrixA=NULL, *matrixB=NULL;	/* Matrices */
static MtxFile_t *fileVin = NULL, 
    *fileVout = NULL;
static uint32_t nocV = 0;               // Tensor product dimension

////////////////////////////////////////////////////////////////////////////////////////////////////

// Convert vector to matrix

static Matrix_t *VecToMat(PTR vec, int fl, int nor, int noc)

{
    Matrix_t *mat;
    PTR d;
    int i, k;

    mat = matAlloc(fl,nor,noc);
    d = mat->data;
    for (i = 0; i < nor; i++) 
    {
        for (k = 0; k < noc; k++) 
            ffInsert(d,k,ffExtract(vec,i*noc+k));
        ffStepPtr(&d, noc);
    }
    return mat;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Convert matrix to vector

static void matToVec(Matrix_t *mat, PTR vec)
{
    const uint32_t noc = mat->noc;
    const uint32_t nor = mat->nor;

    // Clear the vector 
    ffMulRow(vec,FF_ZERO, nor * noc);

    // Convert
    PTR y = mat->data;
    uint32_t l = 0;
    for (uint32_t i = 0; i < nor; ++i) 
    {
        for (uint32_t k = 0; k < noc; ++k) 
            ffInsert(vec,l++,ffExtract(y,k));
        ffStepPtr(&y, noc);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void readMatrices()
{
    matrixA = matLoad(fileNameA);
    matrixB = matLoad(fileNameB);
    if (matrixA == NULL || matrixB == NULL)
    {
	mtxAbort(MTX_HERE,"Error reading matrices");
    }
    if (matrixA->nor != matrixA->noc)
    {
	mtxAbort(MTX_HERE,"%s: %s",fileNameA,MTX_ERR_NOTSQUARE);
    }
    if (matrixB->nor != matrixB->noc)
    {
	mtxAbort(MTX_HERE,"%s: %s",fileNameB,MTX_ERR_NOTSQUARE);
    }
    if (matrixA->field != matrixB->field)
    {
	mtxAbort(MTX_HERE,"%s and %s: %s",fileNameA,fileNameB,MTX_ERR_INCOMPAT);
    }
    nocV = matrixA->noc * matrixB->noc;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void openVectorFiles()
{
    // Open the input vectors file.
    fileVin = mfOpen(fileNameVin);
    mfReadHeader(fileVin);
    if (mfObjectType(fileVin) != MTX_TYPE_MATRIX)
	mtxAbort(MTX_HERE,"%s: %s",fileNameVin,MTX_ERR_NOTMATRIX);
    if (fileVin->header[0] != matrixA->field || fileVin->header[2] != nocV)
	mtxAbort(MTX_HERE,"%s and %s/%s: %s",fileNameVin,fileNameA,fileNameB,MTX_ERR_INCOMPAT);
    
    // Create output file.
    fileVout = mfCreate(fileNameVout,fileVin->header[0],fileVin->header[1],nocV);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char **argv)
{
    App = appAlloc(&AppInfo,argc,argv);
    appGetArguments(App,4,4);
    fileNameVin = App->argV[0];
    fileNameA = App->argV[1];
    fileNameB = App->argV[2];
    fileNameVout = App->argV[3];
    readMatrices();
    openVectorFiles();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void Cleanup()
{
   appFree(App);
   mfClose(fileVout);
   mfClose(fileVin);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{  
    init(argc,argv);

    // Transpose first matrix 
    Matrix_t *matrixATr = matTransposed(matrixA);
    matFree(matrixA);

    // Allocate buffer for one vector
    PTR tmp = ffAlloc(1, nocV);

    // Process input vectors one by one.
    for (uint32_t i = 0; i < fileVin->header[1]; ++i)
    {
	Matrix_t *mat3, *newmat;

        // Read one vector and convert to matrix.
        mfReadRows(fileVin,tmp,1, nocV);
        mat3 = VecToMat(tmp, ffOrder, matrixATr->nor, matrixB->nor);
       
        // Multiply from both sides.
	newmat = matDup(matrixATr);
	matMul(newmat,mat3);
	matMul(newmat,matrixB);

        // Turn matrix into vector and write out
        matToVec(newmat,tmp);
        mfWriteRows(fileVout,tmp,1, nocV);

        matFree(mat3);
        matFree(newmat);
    }
    sysFree(tmp);
   
    Cleanup();
    return 0;
}




////////////////////////////////////////////////////////////////////////////////////////////////////

// *INDENT-OFF*

/**
@page prog_ztm ztm - Tensor Multiply

@section ztm_syntax Command Line
<pre>
ztm [@em Options] @em Vectors @em A @em B @em Result
</pre>

@par @em Options
   Standard options, see @ref prog_stdopts.
@par -s @em Start
    Start with the given seed vector number instead of 1.
@par @em A
    Input file: Left factor (square matrix).
@par @em B
    Input file: Right factor (square matrix).
@par @em Vectors
    Input file: Vectors to be multiplied.
@par @em Result
    Output file: Vectors.

@section ztm_inp Input Files
@par @em Vectors
    Vectors to be multiplied, a r×(mn) matrix.
@par @em A
    Left factor, m×m matrix.
@par @em B
    Right factor, n×n matrix.

@section ztm_out Output Files
@par @em Result
    Vectors.


@section ztm_desc Description
This program reads two matrices from A and B, a list of vectors,
and calculates the image of the vectors under A⊗B.
@em A and @em B must be square matrices of dimension m and n, respectively.

The same calculation could be done with the @ref prog_zte "zte" and @ref prog_zmu "zmu"
programs:
```
ztm Vecs A B Result
```
is equivalent to
```
zte A B AB
zmu Vecs AB Result
```
but using @b ztm avoids the memory-consuming calculation of A⊗B.
**/
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
