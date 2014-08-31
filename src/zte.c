/* ============================= C MeatAxe ==================================
   File:        $Id: zte.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Tensor product of two matrices or permutations.
   --------------------------------------------------------------------------
   (C) Copyright 1999 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

MTX_DEFINE_FILE_INFO

static const char *aname, *bname, *cname;	/* File names */
static MtxFile_t *AFile, *BFile;


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

static int tensormatrices()

{
    PTR m1, m2, m3;
    int a_nor = AFile->Nor;
    int a_noc = AFile->Noc;
    int b_nor = BFile->Nor;
    int b_noc = BFile->Noc;
    int c_nor = a_nor * b_nor;
    int c_noc = a_noc * b_noc;
    int ai;
    MtxFile_t *f;


    MESSAGE(1,("Computing matrix tensor product:"));
    MESSAGE(1,(" (%d,%d)*(%d,%d)=(%d,%d)\n", 
	a_nor,a_noc,b_nor,b_noc,c_nor,c_noc));

    /* Allocate buffers.
       ----------------- */
    FfSetField(AFile->Field);
    FfSetNoc(a_noc);
    m1 = FfAlloc(1);
    FfSetNoc(b_noc);
    m2 = FfAlloc(b_nor);
    FfSetNoc(c_noc);
    m3 = FfAlloc(1);
    if (m1 == NULL || m2 == NULL || m3 == NULL)
	return -1;

    /* Read the second matrix (B).
       --------------------------- */
    FfSetNoc(b_noc);
    if (MfReadRows(BFile,m2,b_nor) != b_nor) 
	return -1;

    /* Open the outpt file.
       -------------------- */
    if ((f = MfCreate(cname,FfOrder,c_nor,c_noc)) == NULL)
	return -1;

    /* Calculate the tensor product.
       ----------------------------- */
    for (ai = 0; ai < a_nor; ++ai)
    {
	int bi;
	PTR bp;

	/* Read the next row from A.
	   ------------------------- */
	if (MfReadRows(AFile,m1,1) != 1)
	    break;

	bp = m2;
	for (bi = 0; bi < b_nor; ++bi)
	{
	    int aj, cj;
	    cj = 0;
	    FfSetNoc(c_noc);
	    for (aj = 0; aj < a_noc; ++aj)
	    {
		FEL f = FfExtract(m1,aj);
		int bj;
		for (bj = 0; bj < b_noc; ++bj)
		{
		    FEL g = FfExtract(bp,bj);
		    FfInsert(m3,cj,FfMul(f,g));
		    ++cj;
		}
	    }
	    if (MfWriteRows(f,m3,1) != 1)
	    {
		MTX_ERROR("Error writing output file");
		return -1;
	    }
	    FfSetNoc(b_noc);
	    FfStepPtr(&bp);
	}
    }
    MfClose(f);
    return 0;
}


/* ------------------------------------------------------------------
    tensorperms() - Calculate the tensor product of two permutations.
   ------------------------------------------------------------------ */

static int tensorperms(void)

{
    int a_deg = AFile->Nor;
    int b_deg = BFile->Nor;
    int c_deg = a_deg * b_deg;
    MtxFile_t *f;
    long *a_buf, *b_buf, *c_buf;
    int i;
    
    MESSAGE(1,("Computing permutation tensor product:"));
    MESSAGE(1,(" %d*%d=%d\n", a_deg,b_deg,c_deg));

    /* Read permutations.
       ------------------ */
    a_buf = NALLOC(long,a_deg);
    b_buf = NALLOC(long,b_deg);
    c_buf = NALLOC(long,b_deg);
    if (a_buf == NULL || b_buf == NULL || c_buf == NULL) 
	return -1;
    if (MfReadLong(AFile,a_buf,a_deg) != a_deg ||
        MfReadLong(BFile,b_buf,b_deg) != b_deg)
    {
	MTX_ERROR("Error reading permutation");
	return -1;
    }

    /* Convert from FORMAT to C format, if necessary.
       ---------------------------------------------- */
    Perm_ConvertOld(a_buf,a_deg);
    Perm_ConvertOld(b_buf,b_deg);

    /* Open output file.
       ----------------- */
    if ((f = MfCreate(cname,-1,c_deg,1)) == NULL)
	return -1;

    /* Calculate the tensor product
       ---------------------------- */
    for (i = 0; i < a_deg; ++i)
    {
	int k;
	for (k = 0; k < b_deg; ++k)
	    c_buf[k] = a_buf[i] * b_deg + b_buf[k];
        if (MfWriteLong(f,c_buf,b_deg) != b_deg)
        {
	    MTX_ERROR("Error writing permutation");
	    return -1;
	}
    }
    MfClose(f);
    return 0;
}


static int OpenFiles()

{
    if ((AFile = MfOpen(aname)) == NULL ||
	(BFile = MfOpen(bname)) == NULL)
    {
	MTX_ERROR("Error opening input file");
	return -1;
    }
    if (AFile->Field != -1 && AFile->Field < 2)
    {
	MTX_ERROR2("%s: Invalid file type %d",aname,AFile->Field);
	return -1;
    }
    if (AFile->Field != BFile->Field)
    {
	MTX_ERROR3("%s and %s: %E",aname,bname,MTX_ERR_INCOMPAT);
	return -1;
    }
    return 0;
}



static int Init(int argc, const char **argv)

{
    /* Process command line options and arguments.
       ------------------------------------------- */
    if ((App = AppAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    if (AppGetArguments(App,3,3) != 3)
	return -1;
    aname = App->ArgV[0];
    bname = App->ArgV[1];
    cname = App->ArgV[2];

    /* Open the input file.
       -------------------- */
    if (OpenFiles() != 0)
	return -1;

    return 0;
}






static void Cleanup()

{
    if (App != NULL)
	AppFree(App);
    if (AFile != NULL)
	MfClose(AFile);
    if (BFile != NULL)
	MfClose(BFile);
}




/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, const char **argv)

{
    int rc = 0;

    if (Init(argc,argv) != 0)
    {
	MTX_ERROR("Initialization failed");
	return 1;
    }
    if (AFile->Field == -1)
	rc = tensorperms();
    else
	rc = tensormatrices();
    Cleanup();
    if (rc != 0) rc = 1;
    return rc;
}

/**
@page prog_zte zte - Tensor Product

<<<<<<< HEAD
@section syntax Command Line
=======
@section zte_syntax Command Line
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
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

<<<<<<< HEAD
@section inp Input Files
=======
@section zte_inp Input Files
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
@par @em A
  Left factor (patrix or permutation).
@par @em B
  Right factor (patrix or permutation).

<<<<<<< HEAD
@section out Output Files
@par @em Result
  Result matrix

@section desc Description
=======
@section zte_out Output Files
@par @em Result
  Result matrix

@section zte_desc Description
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
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
