////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Sum and intersection of two subspaces.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include <stdlib.h>
#include <string.h>
#include "meataxe.h"



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */


static int NorA, NorB;
static int Noc;
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
    if ((of = mfCreate(sumname,ffOrder,NorA,ffNoc)) == NULL)
	return -1;
    if (mfWriteRows(of,Wrk1,NorA) != NorA)
	return -1;
    mfClose(of);

    MESSAGE(1,("Writing intersection to %s\n",intname));
    if ((of = mfCreate(intname,ffOrder,NorB,ffNoc)) == NULL)
	return -1;
    if (mfWriteRows(of,ffGetPtr(Wrk2,NorA,Noc),NorB) != NorB)
	return -1;
    mfClose(of);

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
    if ((af = mfOpen(aname)) == NULL || (bf = mfOpen(bname)) == NULL)
    {
	mtxAbort(MTX_HERE,"Cannot open input file");
	return -1;
    }
    if (af->Field < 2)
    {
	mtxAbort(MTX_HERE,"%s: %s",aname,MTX_ERR_NOTMATRIX);
	return -1;
    }
    if (af->Field != bf->Field || af->Noc != bf->Noc)
    {
	mtxAbort(MTX_HERE,"%s and %s: %s",aname,bname,MTX_ERR_INCOMPAT);
	return -1;
    }
    
    /* Allocate work space.
       -------------------- */
    ffSetField(af->Field);
    Noc = af->Noc;
    ffSetNoc(Noc);
    NorA = af->Nor;
    NorB = bf->Nor;
    Wrk1 = ffAlloc(NorA + NorB, Noc);
    Wrk2 = ffAlloc(NorA + NorB, Noc);
    Piv = NALLOC(int,NorA+NorB);
    if (Wrk1 == NULL || Wrk2 == NULL || Piv == NULL)
	return -1;

    /* Read files.
       ----------- */
    MESSAGE(1,("Reading input files\n"));
    if (mfReadRows(af,Wrk1,NorA) != NorA ||
	mfReadRows(bf,ffGetPtr(Wrk1,NorA,Noc),NorB) != NorB)
    {
	mtxAbort(MTX_HERE,"Error reading input file");
	return -1;
    }

    mfClose(af);
    mfClose(bf);
    return 0;
}



static int Init(int argc, char **argv)

{
    /* Parse command line
       ------------------ */
    if ((App = appAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    if (appGetArguments(App,4,4) < 0)
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
    appFree(App);
}

/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char **argv)
{

    if (Init(argc,argv) != 0)
    {
	mtxAbort(MTX_HERE,"Initialization failed");
	return 1;
    }
    if (ffSumAndIntersection(Noc, Wrk1,&NorA,&NorB,Wrk2,Piv) != 0)
    {
	mtxAbort(MTX_HERE,"Internal error");
	return 1;
    }
    if (WriteFiles() != 0)
	mtxAbort(MTX_HERE,"Error writing output files");
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
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
