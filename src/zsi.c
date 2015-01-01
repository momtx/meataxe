////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Sum and intersection of two subspaces.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include <stdlib.h>
#include <string.h>
#include "meataxe.h"



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

MTX_DEFINE_FILE_INFO

static int NorA, NorB;
static int *Piv;
static PTR Wrk1, Wrk2;
static const char *aname, *bname, *sumname, *intname;


static MtxApplicationInfo_t AppInfo = { 
"zsi", "Sum And Intersection",
"SYNTAX\n"
"    zsi [-QV] <Space1> <Space2> <Sum> <Int>"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"\n"
"ARGUMENTS\n"
"    <Space1> ................ First space\n"
"    <Space2> ................ Second space\n"
"    <Sum> ................... File name for the sum\n"
"    <Int> ................... File name for the intersection\n"
};


static MtxApplication_t *App = NULL;


/* ------------------------------------------------------------------
   WritefFles() - Write out the result
   ------------------------------------------------------------------ */

static int WriteFiles()

{	
    MtxFile_t *of;
    
    MESSAGE(0,("Sum %d, Intersection %d\n",NorA,NorB));
    MESSAGE(1,("Writing sum to %s\n",sumname));
    if ((of = MfCreate(sumname,FfOrder,NorA,FfNoc)) == NULL)
	return -1;
    if (MfWriteRows(of,Wrk1,NorA) != NorA)
	return -1;
    MfClose(of);

    MESSAGE(1,("Writing intersection to %s\n",intname));
    if ((of = MfCreate(intname,FfOrder,NorB,FfNoc)) == NULL)
	return -1;
    if (MfWriteRows(of,FfGetPtr(Wrk2,NorA),NorB) != NorB)
	return -1;
    MfClose(of);

    return 0;
}



/* ------------------------------------------------------------------
   ReadFile() - Read the two spaces.
   ------------------------------------------------------------------ */

static int ReadFiles()

{
    MtxFile_t *af, *bf;

    /* Open files and check headers.
       ----------------------------- */
    if ((af = MfOpen(aname)) == NULL || (bf = MfOpen(bname)) == NULL)
    {
	MTX_ERROR("Cannot open input file");
	return -1;
    }
    if (af->Field < 2)
    {
	MTX_ERROR2("%s: %E",aname,MTX_ERR_NOTMATRIX);
	return -1;
    }
    if (af->Field != bf->Field || af->Noc != bf->Noc)
    {
	MTX_ERROR3("%s and %s: %E",aname,bname,MTX_ERR_INCOMPAT);
	return -1;
    }
    
    /* Allocate work space.
       -------------------- */
    FfSetField(af->Field);
    FfSetNoc(af->Noc);
    NorA = af->Nor;
    NorB = bf->Nor;
    Wrk1 = FfAlloc(NorA + NorB);
    Wrk2 = FfAlloc(NorA + NorB);
    Piv = NALLOC(int,NorA+NorB);
    if (Wrk1 == NULL || Wrk2 == NULL || Piv == NULL)
	return -1;

    /* Read files.
       ----------- */
    MESSAGE(1,("Reading input files\n"));
    if (MfReadRows(af,Wrk1,NorA) != NorA ||
	MfReadRows(bf,FfGetPtr(Wrk1,NorA),NorB) != NorB)
    {
	MTX_ERROR("Error reading input file");
	return -1;
    }

    MfClose(af);
    MfClose(bf);
    return 0;
}



static int Init(int argc, const char **argv)

{
    /* Parse command line
       ------------------ */
    if ((App = AppAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    if (AppGetArguments(App,4,4) < 0)
	return -1;
    aname = App->ArgV[0];
    bname = App->ArgV[1];
    sumname = App->ArgV[2];
    intname = App->ArgV[3];

    if (ReadFiles() != 0)
	return -1;

    return 0;
}


static void Cleanup()

{
    AppFree(App);
}

/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, const char **argv)

{

    if (Init(argc,argv) != 0)
    {
	MTX_ERROR("Initialization failed");
	return 1;
    }
    if (FfSumAndIntersection(Wrk1,&NorA,&NorB,Wrk2,Piv) != 0)
    {
	MTX_ERROR("Internal error");
	return 1;
    }
    if (WriteFiles() != 0)
	MTX_ERROR("Error writing output files");
    Cleanup();
    return 0;
}



/**
@page prog_zsi zsi - Sum and Intersection

@section zsi_syntax Command Line
<pre>
zsi [@em Options] @em Space1 @em Space2 @em Sum @em Int
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par @em Space1
  First space.
@par @em Space2
  Second space.
@par @em Sum
  Sum.
@par @em Int
  Intersection.

@section zsi_inp Input Files
@par @em Space1
  First space.
@par @em Space2
  Second space.

@section zsi_out Output Files
@par @em Sum
  Sum.
@par @em Int
  Intersection.

@section zsi_desc Description

This program reads in two spaces from @em Space1 and @em Space2
and writes out their sum and intersection, in semi-echelon form, to
@em Sum and @em Int, respectively. 
The input files must be matrices over the same field and with the
same number of columns. They need not be in echelon form.

There must be enough memory to hold two copies of each of the two 
spaces at the same time. 
**/
