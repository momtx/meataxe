/* ============================= C MeatAxe ==================================
   File:        $Id: ztc.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Trace of a matrix or permutation.
   --------------------------------------------------------------------------
   (C) Copyright 1999 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include <stdlib.h>



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

MTX_DEFINE_FILE_INFO

static int opt_G = 0;		/* GAP output */
static MtxFile_t *InputFile = NULL;
static const char *inpname = NULL;

static MtxApplicationInfo_t AppInfo = { 
"ztc", "Trace", 
"SYNTAX\n"
"    ztc [-GQV] <File>"
"\n"
"ARGUMENTS\n"
"    <File> .................. Input file name\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -G ...................... GAP output (implies -Q)\n"
"\n"
"FILES\n"
"    <File> .................. I The matrix or permutation\n"
};

static MtxApplication_t *App = NULL;




/* ------------------------------------------------------------------
   trmat() - Trace of a matrix
   ------------------------------------------------------------------ */

static int trmat()

{
    FEL tr;
    int i, max;
    PTR m1;

    FfSetField(InputFile->Field); 
    FfSetNoc(InputFile->Noc);
    m1 = FfAlloc(1);
    tr = FF_ZERO;
    if ((max = InputFile->Nor) > InputFile->Noc) 
	max = InputFile->Noc;
    for (i = 0; i < max; ++i)
    {
	if (MfReadRows(InputFile,m1,1) != 1)
	{
	    MTX_ERROR("Cannot read ipnut file");
	    return -1;
	}
	tr = FfAdd(tr,FfExtract(m1,i));
    }
    if (!opt_G)		/* Standard output */
	printf("Trace is %d\n",FfToInt(tr));
    else		/* GAP output */
        printf("MeatAxe.Trace := %s;\n",FfToGap(tr));
    return 0;
}


/* ------------------------------------------------------------------
   trperm() - Trace of a permutation (= number of fixed points)
   ------------------------------------------------------------------ */

static int trperm()

{
    long *m1, tr;
    int degree = InputFile->Nor;
    int k;

    m1 = NALLOC(long,degree);
    if (m1 == NULL) 
	return -1;

    if (MfReadLong(InputFile,m1,degree) != degree)
    {
	MTX_ERROR("Error reading permutation");
	return -1;
    }

    tr = 0;
    for (k = 0; k < degree; ++k)
    {
	if (m1[k] == k) 
	    ++tr;
    }

    if (!opt_G)
	printf("Trace is %ld\n",tr);
    else
	printf("MeatAxe.Trace := [%ld];\n",tr);
    return 0;
}


static int Init(int argc, const char **argv)

{
    /* Process command line options and arguments.
       ------------------------------------------- */
    if ((App = AppAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    opt_G = AppGetOption(App,"-G --gap");
    if (opt_G) MtxMessageLevel = -100;
    if (AppGetArguments(App,1,1) != 1)
	return -1;
    inpname = App->ArgV[0];

    /* Open the input file.
       -------------------- */
    if ((InputFile = MfOpen(inpname)) == NULL)
    {
	MTX_ERROR("Error opening input file");
	return -1;
    }
    if (InputFile->Field >= 2)
    {
        MESSAGE(1,("Input is a %dx%d matrix over GF(%d)\n",InputFile->Nor,
	    InputFile->Noc,InputFile->Field));
    }
    else if (InputFile->Field == -1)
    {
        MESSAGE(1,("Input is a permutation of degree %d\n",InputFile->Nor));
    }
    else
    {
	MTX_ERROR2("%s: Unknown type %d",inpname,InputFile->Field);
	return -1;
    }

    return 0;
}


static void Cleanup()

{
    if (App != NULL)
	AppFree(App);
    if (InputFile != NULL)
	MfClose(InputFile);
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
    if (InputFile->Field == -1) 
	rc = trperm() != 0 ? 1 : 0;
    else
	rc = trmat() != 0 ? 1 : 0;
    Cleanup();
    return rc;
}



/**
@page prog_ztc ztc - Trace

<<<<<<< HEAD
@section syntax Command Line
=======
@section ztc_syntax Command Line
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
<pre>
ztc [@em Options] [-G] @em Inp
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par -G
  GAP output.
@par @em Inp
  Input matrix or permutation.

<<<<<<< HEAD
@section inp Input Files
@par @em Inp
  Input matrix or permutation.

@section desc Description
=======
@section ztc_inp Input Files
@par @em Inp
  Input matrix or permutation.

@section ztc_desc Description
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
This program reads a matrix or permutation calculates its trace and
outputs this to the user.
If the input file is a matrix, it is read row by row and the diagonal
entries are added up. Note that the matrix need not be square. The 
result is printed in the form
<pre>
Trace is N
</pre>
For permutations, the trace is calculated as the number of fixed points.
*/
