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
"    ztm [-QV] [-T <MaxTime>] <A> <B> <Vectors> <Result>"
"\n"
"ARGUMENTS\n"
"    <A> ..................... Input file: Left factor (square matrix)\n"
"    <B> ..................... Input file: Right factor (square matrix)\n"
"    <Vectors> ............... Input file: Vectors\n"
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



/* ------------------------------------------------------------------
   VecToMat() - Convert vector to matrix
   ------------------------------------------------------------------ */

static Matrix_t *VecToMat(PTR vec, int fl, int nor, int noc)

{
    Matrix_t *mat;
    PTR d;
    int i, k;

    mat = matAlloc(fl,nor,noc);
    d = mat->Data;
    for (i = 0; i < nor; i++) 
    {
        for (k = 0; k < noc; k++) 
            ffInsert(d,k,ffExtract(vec,i*noc+k));
        ffStepPtr(&d, noc);
    }
    return mat;
}




/* ------------------------------------------------------------------
   matToVec() - Convert matrix to vector
   ------------------------------------------------------------------ */

static void matToVec(Matrix_t *mat, PTR vec)

{
    PTR y;
    int i, k, nrows, ncol, l;

    ncol = mat->Noc;
    nrows = mat->Nor;

    /* Clear the vector 
       ---------------- */
    ffMulRow(vec,FF_ZERO, nrows * ncol);

    /* Convert
       ------- */
    y = mat->Data;
    l = 0;
    for (i = 0; i < nrows; ++i) 
    {
        for (k = 0; k < ncol; ++k) 
            ffInsert(vec,l++,ffExtract(y,k));
        ffStepPtr(&y, ncol);
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
    if (matrixA->Nor != matrixA->Noc)
    {
	mtxAbort(MTX_HERE,"%s: %s",fileNameA,MTX_ERR_NOTSQUARE);
    }
    if (matrixB->Nor != matrixB->Noc)
    {
	mtxAbort(MTX_HERE,"%s: %s",fileNameB,MTX_ERR_NOTSQUARE);
    }
    if (matrixA->Field != matrixB->Field)
    {
	mtxAbort(MTX_HERE,"%s and %s: %s",fileNameA,fileNameB,MTX_ERR_INCOMPAT);
    }
    nocV = matrixA->Noc * matrixB->Noc;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void openVectorFiles()
{
    // Open the input vectors file.
    fileVin = mfOpen(fileNameVin);
    mfReadHeader(fileVin);
    if (mfObjectType(fileVin) != MTX_TYPE_MATRIX)
	mtxAbort(MTX_HERE,"%s: %s",fileNameVin,MTX_ERR_NOTMATRIX);
    if (fileVin->header[0] != matrixA->Field || fileVin->header[2] != nocV)
	mtxAbort(MTX_HERE,"%s and %s/%s: %s",fileNameVin,fileNameA,fileNameB,MTX_ERR_INCOMPAT);
    
    // Create output file.
    fileVout = mfCreate(fileNameVout,fileVin->header[0],fileVin->header[1],nocV);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char **argv)
{
    App = appAlloc(&AppInfo,argc,argv);
    appGetArguments(App,4,4);
    fileNameA = App->ArgV[0];
    fileNameB = App->ArgV[1];
    fileNameVin = App->ArgV[2];
    fileNameVout = App->ArgV[3];
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
        mat3 = VecToMat(tmp, ffOrder, matrixATr->Nor, matrixB->Nor);
       
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




/**
@page prog_ztm ztm - Tensor Multiply

@section ztm_syntax Command Line
<pre>
ztm [@em Options] @em A @em B @em Vectors @em Result
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
@par @em A
    Left factor, m×m matrix.
@par @em B
    Right factor, n×n matrix.
@par @em Vectors
    Vectors to be multiplied, a r×(mn) matrix.

@section ztm_out Output Files
@par @em Result
    Vectors.


@section ztm_desc Description
This program reads two matrices from A and B, a list of vectors,
and calculates the image of the vectors under A⊗B.
@em A and @em B must be square matrices of dimension m and n, respetively
This calculation could be done with the @ref prog_zte "zte" and @ref prog_zmu "zmu"
programs, but using @b ztm avoids the memory-consuming calculation of A⊗B.

@subsection ztm_impl Implementation Details
Let m, n be the dimensions of A and B, respectively. For each
input vector v∊F<sup>1×mn</sup>, the program computes the m×n matrix \f$\hat{v}\f$ with
@f[
	\hat{v}_{ik} = v_{(i-1)n + k}\qquad (1\leq i\leq m,\ 1\leq k\leq n)
@f]
Then, v(A⊗B) can be computed by ordinary matrix operations:
@f[
    \widehat{v(A\otimes B)} = A^{tr}\cdot\hat{v}\cdot B
@f]

Both matrices and one vector must fit into memory at the same time.


**/
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
