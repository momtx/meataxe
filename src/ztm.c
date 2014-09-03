/* ============================= C MeatAxe ==================================
   File:        $Id: ztm.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Map vector under tensor product.
   --------------------------------------------------------------------------
   (C) Copyright 1997 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include <stdlib.h>
#include <string.h>


/* ------------------------------------------------------------------
   Variables
   ------------------------------------------------------------------ */

MTX_DEFINE_FILE_INFO

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
static const char *AName, *BName;	/* Matrices */
static const char *VName, *RName;	/* Vectors, Result */
static Matrix_t *mat1=NULL, *mat2=NULL;	/* Matrices */
static MtxFile_t *vecfile = NULL, 
    *resultfile = NULL;


/* ------------------------------------------------------------------
   VecToMat() - Convert vector to matrix
   ------------------------------------------------------------------ */

static Matrix_t *VecToMat(PTR vec, int fl, int nor, int noc)

{
    Matrix_t *mat;
    PTR d;
    int i, k;

    mat = MatAlloc(fl,nor,noc);
    d = mat->Data;
    FfSetNoc(noc);
    for (i = 0; i < nor; i++) 
    {
        for (k = 0; k < noc; k++) 
            FfInsert(d,k,FfExtract(vec,i*noc+k));
        FfStepPtr(&d);
    }
    return mat;
}




/* ------------------------------------------------------------------
   MatToVec() - Convert matrix to vector
   ------------------------------------------------------------------ */

static void MatToVec(Matrix_t *mat, PTR vec)

{
    PTR y;
    int i, k, nrows, ncol, l;

    ncol = mat->Noc;
    nrows = mat->Nor;

    /* Clear the vector 
       ---------------- */
    FfSetNoc(ncol * nrows);
    FfMulRow(vec,FF_ZERO);

    /* Convert
       ------- */
    FfSetNoc(ncol);
    y = mat->Data;
    l = 0;
    for (i = 0; i < nrows; ++i) 
    {
        for (k = 0; k < ncol; ++k) 
            FfInsert(vec,l++,FfExtract(y,k));
        FfStepPtr(&y);
    }
}




static int ReadMatrices()
{
    /* Read matrices and check compatibility.
       -------------------------------------- */
    mat1 = MatLoad(AName);
    mat2 = MatLoad(BName);
    if (mat1 == NULL || mat2 == NULL)
    {
	MTX_ERROR("Error reading matrices");
	return -1;
    }
    if (mat1->Nor != mat1->Noc)
    {
	MTX_ERROR2("%s: %E",AName,MTX_ERR_NOTSQUARE);
	return -1;
    }
    if (mat2->Nor != mat2->Noc)
    {
	MTX_ERROR2("%s: %E",BName,MTX_ERR_NOTSQUARE);
	return -1;
    }
    if (mat1->Field != mat2->Field)
    {
	MTX_ERROR3("%s and %s: %E",AName,BName,MTX_ERR_INCOMPAT);
	return -1;
    }
    MESSAGE(1,("%s: %dx%d matrix over GF(%d)\n", AName,mat1->Nor,mat1->Noc,mat1->Field));
    MESSAGE(1,("%s: %dx%d matrix over GF(%d)\n", BName,mat2->Nor,mat2->Noc,mat2->Field));
    return 0;
}


static int OpenVectorFiles()
{
    /* Open the <Vectors> file and check the header.
       --------------------------------------------- */
    vecfile = MfOpen(VName);
    if (vecfile == NULL)
	return -1;
    if (vecfile->Field != mat1->Field 
	|| vecfile->Noc != mat1->Noc * mat2->Noc)
    {
	MTX_ERROR4("%s and %s/%s: %E",VName,AName,BName,MTX_ERR_INCOMPAT);
	return -1;
    }
    MESSAGE(1,("%s: %dx%d matrix over GF(%d)\n",VName,vecfile->Nor,vecfile->Noc,vecfile->Field));
    
    /* Create output file.
       ------------------- */
    resultfile = MfCreate(RName,vecfile->Field,vecfile->Nor,vecfile->Noc);
    if (resultfile == NULL)
	return -1;

    return 0;
}


static int Init(int argc, const char **argv)

{

    /* Process command line options and arguments.
       ------------------------------------------- */
    if ((App = AppAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    if (AppGetArguments(App,4,4) != 4)
	return -1;
    AName = App->ArgV[0];
    BName = App->ArgV[1];
    VName = App->ArgV[2];
    RName = App->ArgV[3];

    if (ReadMatrices() != 0)
	return -1;
    if (OpenVectorFiles() != 0)
	return -1;

    return 0;
}



static void Cleanup()

{
    if (App != NULL)
	AppFree(App);
    if (resultfile != NULL)
	MfClose(resultfile);
    if (vecfile != NULL)
	MfClose(vecfile);
}



/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, const char **argv)

{  
    Matrix_t *mat1tr;
    PTR tmp;
    long i;


    if (Init(argc,argv) != 0)
    {
	MTX_ERROR("Initialization failed");
	return 1;
    }

    /* Transpose first matrix 
       ---------------------- */
    mat1tr = MatTransposed(mat1);
    MatFree(mat1);

    /* Allocate buffer for one vector
       ------------------------------ */
    FfSetNoc(vecfile->Noc);
    tmp = FfAlloc(1);

    /* Process <Vectors> one by one
       ---------------------------- */
    for (i = 0; i < vecfile->Nor; ++i)
    {
	Matrix_t *mat3, *newmat;

        /* Read on evector and convert to matrix.
           -------------------------------------- */
        if (MfReadRows(vecfile,tmp,1) != 1)
	{
	    MTX_ERROR("Error reading vector");
	    return 1;
	}
        mat3 = VecToMat(tmp,FfOrder,mat1tr->Nor,mat2->Nor);
       
        /* Multiply from both sides.
           ------------------------- */
	newmat = MatDup(mat1tr);
	MatMul(newmat,mat3);
	MatMul(newmat,mat2);

        /* Turn matrix into vector and write out
           ------------------------------------- */
        MatToVec(newmat,tmp);
        if (MfWriteRows(resultfile,tmp,1) != 1)
	{
	    MTX_ERROR("Error writing vector");
	    return 1;
	}

        /* Free memory
           ----------- */
        MatFree(mat3);
        MatFree(newmat);
    }
   
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
