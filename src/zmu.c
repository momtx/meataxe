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
"    zmu [-r <Row>[.<#Rows>]] [-c <Col>[.<#Cols>]] <A> <B> <Result>\n"
"\n"
"OPTIONS\n"
"    -r and -c can be used for piecewise matrix multiplication\n"
"    E.g., `-r1 -c2' selects the upper right part. By default,\n"
"    <#Rows> = <#Cols> = 2, i.e., the result is divided into\n"
"    four pieces.\n"
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


static MtxApplication_t *App = NULL;
static const char *aname, *bname, *cname;	/* File names */
static FILE *afile = NULL; 
static FILE *bfile = NULL;
static FILE *cfile = NULL;		    /* File handles */
static int fl1, fl2;			    /* Field order */
static int nor1, noc1, nor2, noc2;	    /* Matrix dimensions */
int nrows = 1, thisrow = 1;		    /* Arguments to -r option */
int ncols = 1, thiscol = 1;		    /* Arguments to -c option */




/* ------------------------------------------------------------------
   multpm() - Multiply permutation by matrix
   ------------------------------------------------------------------ */

static int multpm(void)

{
    int i;
    Perm_t *perm;
    Matrix_t *row;

    if (nor1 != nor2) 
    {
	mtxAbort(MTX_HERE,"%s and %s: %s",aname,bname,MTX_ERR_INCOMPAT);
	return -1;
    }

    /* Read the permutation (A).
       ------------------------- */
    sysFseek(afile,0);
    perm = permRead(afile);
    if (perm == NULL)
    {
	mtxAbort(MTX_HERE,"Cannot read permutation from %s",aname);
	return -1;
    }

    /* Allocate workspace (one row of <B>).
       ------------------------------------ */
    row = matAlloc(fl2,1,noc2);

    /* Open the output file.
       --------------------- */
    if ((cfile = ffWriteHeader(cname,fl2,nor2,noc2)) == NULL)
	return -1;

    /* Write out the rows of <B> in the order defined by <A>.
       ------------------------------------------------------ */
    for (i = 0; i < nor1; ++i)
    {	
	ffSeekRow(bfile,perm->Data[i], noc2);
	ffReadRows(bfile,row->Data, 1, noc2);
	ffWriteRows(cfile,row->Data,1, noc2);
    }

    /* Clean up.
       --------- */
    matFree(row);
    permFree(perm);
    return 0;
}


/* ------------------------------------------------------------------
   multmp() - Multiply matrix by permutation
   ------------------------------------------------------------------ */

static int multmp(void)

{
    Perm_t *perm;
    Matrix_t *buffer;
    PTR row_in, row_out;
    int i;

    if (noc1 != nor2) 
    {
	mtxAbort(MTX_HERE,"%s and %s: %s",aname,bname,MTX_ERR_INCOMPAT);
	return -1;
    }

    /* Read the permutation (B).
       ------------------------- */
    sysFseek(bfile,0);
    perm = permRead(bfile);
    if (perm == NULL)
    {
	mtxAbort(MTX_HERE,"Cannot read permutation from %s",bname);
	return -1;
    }

    /* Allocate workspace (two rows of A).
       ----------------------------------- */
    buffer = matAlloc(fl1,2,noc1);
    row_in = matGetPtr(buffer,0);
    row_out = matGetPtr(buffer,1);

    /* Create the output file.
       ----------------------- */
    if ((cfile = ffWriteHeader(cname,fl1,nor1,noc1)) == NULL)
	return -1;

    /* Process A row by row. Permute the 
       marks of each row according to B.
       ---------------------------------- */
    for (i = 0; i < nor1; ++i)
    {
        ffReadRows(afile,row_in,1, noc1);
	ffPermRow(row_in,perm->Data,noc1,row_out);
	ffWriteRows(cfile,row_out,1, noc1);
    }

    /* Clean up.
       --------- */
    matFree(buffer);
    permFree(perm);
    return 0;
}



/* ------------------------------------------------------------------
   multsm() - Multiply scalar with matrix
   ------------------------------------------------------------------ */

static int multsm(FILE *s, FILE *m, int norM, int nocM)
{
    ffSetField(fl1);
    PTR ms = ffAlloc(1,1);
    ffReadRows(s,ms,1,1);
    const FEL f = ffExtract(ms,0);
    PTR mm = ffAlloc(1, nocM);
	
    if ((cfile = ffWriteHeader(cname,fl1,norM,nocM)) == NULL)
	return -1;
    for (int i = 0; i < norM; ++i)
    {	
	ffReadRows(m,mm,1, nocM);
	ffMulRow(mm,f, nocM);
	ffWriteRows(cfile,mm,1, nocM);
    }
    sysFree(mm);
    return 0;
}



/* ------------------------------------------------------------------
   multmm() - Multiply two matrices
   ------------------------------------------------------------------ */

static int multmm(void)
{
    PTR m1, m2, tmp, out;
    int row1, row2;	/* Range of rows */
    int col1, col2;	/* Range of columns */
    int i, k;
    int diffsize;

    if (fl1 != fl2)
	mtxAbort(MTX_HERE,"%s and %s: %s (different fields)",aname,bname, MTX_ERR_INCOMPAT);

    if (noc1 == 1 && nor1 == 1)
	return multsm(afile,bfile,nor2,noc2);
    else if (noc2 == 1 && nor2 == 1)
	return multsm(bfile,afile,nor1,noc1);
    if (noc1 != nor2) 
	mtxAbort(MTX_HERE,"%s and %s: %s",aname,bname,MTX_ERR_INCOMPAT);
    if ((long) nrows > nor1 || (long) ncols > noc2) 
	mtxAbort(MTX_HERE,"Matrix too small");
    row1 = (nor1 / nrows) * (thisrow - 1);
    row2 = (nor1 / nrows) * thisrow - 1;
    col1 = (noc2 / ncols) * (thiscol - 1);
    col2 = (noc2 / ncols) * thiscol - 1;

    /* First matrix
       ------------ */
    ffSetField(fl1);
    m1 = ffAlloc(1, noc1);
    ffSeekRow(afile,row1, noc1);

    /* Read second matrix
       ------------------ */
    tmp = ffAlloc(1, noc2);
    const int colspan = col2 - col1 + 1;
    diffsize = ffRowSize(colspan);
    m2 = ffAlloc(nor2, colspan);
    out = ffAlloc(1, colspan);
    if (colspan < noc2)
    {
	PTR x = m2;
	for (i = 1; i <= nor2; ++i)
	{
	    ffReadRows(bfile, tmp, 1, noc2);
	    for (k = col1; k <= col2; ++k)
		ffInsert(x,k - col1,ffExtract(tmp,k));
	    x = (PTR)((char *)x + diffsize);
	}
    }
    else
    {
	ffReadRows(bfile,m2,nor2, noc2);
    }

    /* Multiply and write result
       ------------------------- */
    if ((cfile = ffWriteHeader(cname,fl1,row2-row1+1,col2-col1+1)) == NULL)
	return -1;

    const int nocOut = col2 - col1 + 1;
    for (i = 0; i < row2 - row1 + 1; ++i)
    {
        ffReadRows(afile, m1, 1, noc1);
	ffMapRow(m1, m2, nor2, nocOut, out);
	ffWriteRows(cfile,out,1, nocOut);
    }

    return 0;
}


/* ------------------------------------------------------------------
   multpp() - Multiply two permutations
   ------------------------------------------------------------------ */

static int multpp(void)

{
    Perm_t *a, *b;

    if (fl1 != fl2) 
    {
	mtxAbort(MTX_HERE,"%s and %s: %s",aname,bname,MTX_ERR_INCOMPAT);
	return -1;
    }
    if (fl1 != -1) 
    {
	mtxAbort(MTX_HERE,"Monomials are not supported");
	return -1;
    }

    /* Read the permutations.
       ---------------------- */
    sysFseek(afile,0);
    a = permRead(afile);
    sysFseek(bfile,0);
    b = permRead(bfile);

    /* Multiply and write.
       ------------------- */
    permMul(a,b);
    permSave(a,cname);

    /* Clean up.
       --------- */
    permFree(a);
    permFree(b);
    return 0;    
}



static int ParseSpec(const char *spec, const char *name, int *pos, int *size)
{
    *pos = atoi(spec);
    while (*spec != 0 && isdigit(*spec)) ++spec;
    if (*spec == '.' || *spec == '/')
	*size = atoi(spec+1);
    else
	*size = 2;

    if (*pos < 1 || *pos > *size)
    {
	mtxAbort(MTX_HERE,"Invalid %s specification",name);
	return -1;
    }
    return 0;
}


static int Init(int argc, char **argv)
{
    const char *c;

    App = appAlloc(&AppInfo,argc,argv);
    if (App == NULL)
	return -1;

    /* Command line options.
       --------------------- */
    if ((c = appGetTextOption(App,"-c",NULL)) != NULL)
	ParseSpec(c,"column",&thiscol,&ncols);
    if ((c = appGetTextOption(App,"-r",NULL)) != NULL)
	ParseSpec(c,"row",&thisrow,&nrows);
    if (appGetArguments(App,3,3) < 0)
	return -1;

    /* Command line arguments.
       ----------------------- */
    aname = App->ArgV[0];
    bname = App->ArgV[1];
    cname = App->ArgV[2];
    if (!strcmp(aname,cname) || !strcmp(bname,cname))
    {
	mtxAbort(MTX_HERE,"Identical file names not allowed");
	return -1;
    }

    return 0;

}


static int OpenFiles()

{
    afile = ffReadHeader(aname,&fl1,&nor1,&noc1);
    if (afile == NULL) 
	return -1;
    bfile = ffReadHeader(bname,&fl2,&nor2,&noc2);
    if (bfile == NULL) 
	return -1;
    return 0;
}



static void Cleanup()

{
    if (afile != NULL)
	fclose(afile);
    if (bfile != NULL)
        fclose(bfile);
    if (cfile != NULL)
	fclose(cfile);
    appFree(App);
}


int main(int argc, char **argv)
{   
    int result;

    if (Init(argc,argv) != 0)
	return -1;

    if (OpenFiles() != 0)
	return -1;

    /* Call the appropriate multiplication function.
       --------------------------------------------- */
    if (fl1 < 0 && fl2  < 0) 
	result = multpp(); 
    else if (fl1 > 1 && fl2 > 1) 
	result = multmm(); 
    else if (fl1 > 1 && fl2 == -1) 
	result = multmp(); 
    else if (fl1 == -1 && fl2  > 1) 
	result = multpm(); 
    else
    {
	mtxAbort(MTX_HERE,"%s and %s: %s",aname,bname,MTX_ERR_INCOMPAT);
	result = 1;
    }

    Cleanup();

    return result;
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
