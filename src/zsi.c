////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Sum and intersection of two subspaces.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include <stdlib.h>
#include <string.h>
#include "meataxe.h"



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */


static uint32_t norA, norB;
static uint32_t noc;
static uint32_t *Piv;
static PTR wrk1, wrk2;
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

////////////////////////////////////////////////////////////////////////////////////////////////////

static void writeFiles()
{	
    MESSAGE(0,("Sum %d, Intersection %d\n",norA,norB));

    MtxFile_t* of = mfCreate(sumname,ffOrder,norA,noc);
    mfWriteRows(of,wrk1,norA, noc);
    mfClose(of);

    of = mfCreate(intname,ffOrder,norB,noc);
    mfWriteRows(of,ffGetPtr(wrk2,norA,noc),norB, noc);
    mfClose(of);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void readFiles()
{
    MtxFile_t* af = mfOpen(aname);
    mfReadHeader(af);
    if (mfObjectType(af) != MTX_TYPE_MATRIX)
       mtxAbort(MTX_HERE, "%s: %s", aname, MTX_ERR_NOTMATRIX);
    norA = af->header[1];
    noc = af->header[2];

    MtxFile_t* bf = mfOpen(bname);
    mfReadHeader(bf);
    if (mfObjectType(bf) != MTX_TYPE_MATRIX)
       mtxAbort(MTX_HERE, "%s: %s", bname, MTX_ERR_NOTMATRIX);
    norB = bf->header[1];
    if (bf->header[0] != af->header[0] || bf->header[2] != noc)
	mtxAbort(MTX_HERE,"%s and %s: %s",aname,bname,MTX_ERR_INCOMPAT);
    
    // Allocate work spaces.
    ffSetField(af->header[0]);
    wrk1 = ffAlloc(norA + norB, noc);
    wrk2 = ffAlloc(norA + norB, noc);
    Piv = NALLOC(uint32_t,norA+norB);

    // Read both subspaces into wrk1.
    mfReadRows(af,wrk1,norA,noc);
    mfReadRows(bf,ffGetPtr(wrk1,norA,noc),norB, noc);

    mfClose(af);
    mfClose(bf);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char **argv)
{
    App = appAlloc(&AppInfo,argc,argv);
    appGetArguments(App,4,4);
    aname = App->argV[0];
    bname = App->argV[1];
    sumname = App->argV[2];
    intname = App->argV[3];
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    init(argc,argv);
    readFiles();
    ffSumAndIntersection(noc, wrk1, &norA, &norB, wrk2, Piv);
    writeFiles();
    appFree(App);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// *INDENT-OFF*

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
