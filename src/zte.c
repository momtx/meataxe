////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tensor product of two matrices or permutations.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"


static const char *fileNameA, *fileNameB, *fileNameC;
static MtxFile_t *fileA, *fileB;


static MtxApplicationInfo_t AppInfo = { 
"zte", "Tensor Product", 
"SYNTAX\n"
"    zte [-QV] [-T <MaxTime>] <A> <B> <Result>"
"\n"
"ARGUMENTS\n"
"    <A> ..................... Left factor\n"
"    <B> ..................... Right factor\n"
"    <Result> ................ Tensor product\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
};

static MtxApplication_t *App = NULL;


/* ------------------------------------------------------------------
   tensormatrices() - Calculate the tensor product of two matrices
   ------------------------------------------------------------------ */

static void tensormatrices()
{
    const uint32_t field = fileA->header[0];
    const uint32_t norA = fileA->header[1];
    const uint32_t nocA = fileA->header[2];
    const uint32_t norB = fileB->header[1];
    const uint32_t nocB = fileB->header[2];
    const uint32_t norC = norA * norB;
    const uint32_t nocC = nocA * nocB;


    MESSAGE(1, "Computing matrix tensor product:");
    MESSAGE(1, " (%d,%d)*(%d,%d)=(%d,%d)\n", norA,nocA,norB,nocB,norC,nocC);

    // Allocate buffers.
    ffSetField(field);
    PTR m1 = ffAlloc(1, nocA);
    PTR m2 = ffAlloc(norB, nocB);
    PTR m3 = ffAlloc(1, nocC);

    // Read the second matrix (B).
    mfReadRows(fileB,m2,norB,nocB);

    // Open the outpt file.
    MtxFile_t *fileC = mfCreate(fileNameC,ffOrder,norC,nocC);

    // Calculate the tensor product.
    for (uint32_t ai = 0; ai < norA; ++ai)
    {
	int bi;
	PTR bp;

	// Read the next row from A.
	mfReadRows(fileA,m1,1,norA);

	bp = m2;
	for (bi = 0; bi < norB; ++bi)
	{
	    int aj, cj;
	    cj = 0;
	    for (aj = 0; aj < nocA; ++aj)
	    {
		FEL f = ffExtract(m1,aj);
		int bj;
		for (bj = 0; bj < nocB; ++bj)
		{
		    FEL g = ffExtract(bp,bj);
		    ffInsert(m3,cj,ffMul(f,g));
		    ++cj;
		}
	    }
	    mfWriteRows(fileC,m3,1,nocC);
	    ffStepPtr(&bp, nocB);
	}
    }
    mfClose(fileC);
}


/* ------------------------------------------------------------------
    tensorperms() - Calculate the tensor product of two permutations.
   ------------------------------------------------------------------ */

static int tensorperms(void)
{
    const uint32_t a_deg = fileA->header[1];
    const uint32_t b_deg = fileB->header[1];
    MTX_ASSERT((((uint64_t)a_deg * b_deg) >> 32) == 0);
    const uint32_t c_deg = a_deg * b_deg;
    int i;
    
    MESSAGE(1, "Computing permutation tensor product:");
    MESSAGE(1, " %d*%d=%d\n", a_deg,b_deg,c_deg);

    // Read permutations.
    uint32_t* a_buf = NALLOC(uint32_t,a_deg);
    uint32_t* b_buf = NALLOC(uint32_t,b_deg);
    uint32_t* c_buf = NALLOC(uint32_t,b_deg);
    sysRead32(fileA->file,a_buf,a_deg);
    sysRead32(fileB->file,b_buf,b_deg);
    permConvertLegacyFormat(a_buf,a_deg);
    permConvertLegacyFormat(b_buf,b_deg);

    // Open output file.
    MtxFile_t *f = mfCreate(fileNameC,-1,c_deg,1);

    // Calculate the tensor product
    for (i = 0; i < a_deg; ++i)
    {
	int k;
	for (k = 0; k < b_deg; ++k)
	    c_buf[k] = a_buf[i] * b_deg + b_buf[k];
        sysWrite32(f->file, c_buf, b_deg);
    }
    mfClose(f);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void doit()
{
    fileA = mfOpen(fileNameA);
    mfReadHeader(fileA);
    const uint32_t objectType = mfObjectType(fileA);
    if (objectType != MTX_TYPE_MATRIX && objectType != MTX_TYPE_PERMUTATION) {
	mtxAbort(MTX_HERE,"%s: unsupported object type 0x%lx",
              fileNameA,(unsigned long) objectType);
    }

    fileB = mfOpen(fileNameB);
    mfReadHeader(fileB);
    if (mfObjectType(fileB) != objectType) {
	mtxAbort(MTX_HERE,"%s and %s: %s",fileNameA,fileNameB,MTX_ERR_INCOMPAT);
    }
    if (objectType == MTX_TYPE_MATRIX) {
       if (fileB->header[0] != fileA->header[0]) {
          mtxAbort(MTX_HERE,"%s and %s: %s (different fields)",
                fileNameA,fileNameB,MTX_ERR_INCOMPAT);
       }
       tensormatrices();
    } else {
       tensorperms();
    }
    mfClose(fileA);
    mfClose(fileB);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char **argv)
{
    App = appAlloc(&AppInfo,argc,argv);
    appGetArguments(App,3,3);
    fileNameA = App->argV[0];
    fileNameB = App->argV[1];
    fileNameC = App->argV[2];
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void cleanup()
{
   appFree(App);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    init(argc,argv);
    doit();
    cleanup();
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// *INDENT-OFF*

/**
@page prog_zte zte - Tensor Product

@section zte_syntax Command Line
<pre>
zte [@em Options] @em A @em B @em Result
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par @em A
  Left factor.
@par @em B
  Right factor.
@par @em Result
  Result matrix

@section zte_inp Input Files
@par @em A
  Left factor (patrix or permutation).
@par @em B
  Right factor (patrix or permutation).

@section zte_out Output Files
@par @em Result
  Result matrix

@section zte_desc Description
This program reads in two matrices or permutations and writes out their tensor
(Kronecker) product.
If @em A is an m×n matrix, and @em B is an m'×n' matrix, the
result is an mm'×nn' matrix consisting of all possible products
of entries from @em A and @em B.
An example may make this clearer:
<pre>
$ zpr G1
matrix field=17 nor=2 noc=3
1 2 3
1 1 0
$ zpr G2
matrix field=17 nor=2 noc=2
2 0
1 3
$ zte G1 G2 P1
$ zpr P1
matrix field=17 nor=2 noc=6
2 0 4 0 6 0
1 3 2 6 3 9
2 0 2 0 0 0
1 3 1 3 0 0
</pre>
If the input files contain permutations on n and n' points, 
respectively, the result is a permutation on nn' points, which 
describes the action of (A,B) on ordered pairs of points. This 
action is defined in the obvious way: (i,k) maps to (iA,kB). 
In the output, pairs are represented as numbers using the 
lexicographic ordering (1,1), ..., (1,n'), (2,1), ..., (n,n').
**/
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
