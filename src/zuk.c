////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Uncondense vectors
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"
#include <stdlib.h>
#include <string.h>


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */



static MtxApplicationInfo_t AppInfo = { 
"zuk", "Uncondense Vectors", 
"SYNTAX\n"
"    zuk " MTX_COMMON_OPTIONS_SYNTAX " <Vectors> <Orbits> <Result>\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"\n"
"FILES\n"
"    <Vectors> ............... I A matrix (#columns = #orbits)\n"
"    <Orbits> ................ I Orbits file, produces by ZMO\n"
"    <Result> ................ O Uncondensed vectors\n"
};

static MtxApplication_t *App = NULL;
static const char *vecname, *orbname, *resname;
static IntMatrix_t *Orbits = NULL;
static IntMatrix_t *OrbitSizes = NULL;
static MtxFile_t *InputFile = NULL;
static MtxFile_t *OutputFile = NULL;
static PTR InputBuffer = NULL;
static PTR OutputBuffer = NULL;
static int Degree, NOrbits;


/* ------------------------------------------------------------------
   uncondense()
   ------------------------------------------------------------------ */

static int uncondense()

{
    long i;

    /* Uncondense row by row
       --------------------- */
    for (i = 0; i < InputFile->Nor; ++i)
    {
	long k;

	if (mfReadRows(InputFile,InputBuffer,1) != 1)
	{
	    mtxAbort(MTX_HERE,"Error reading vector from %s",InputFile->Name);
	    return -1;
	}

	ffMulRow(OutputBuffer,FF_ZERO, Degree);

	for (k = 0; k < NOrbits; ++k)
	{
	    int l;
	    FEL f = ffExtract(InputBuffer,k);
	    int count = OrbitSizes->Data[k];
	    for (l = 0; count > 0 && l < Degree; ++l)
	    {
		if (Orbits->Data[l] == k)
		{
		    ffInsert(OutputBuffer,l,f);
		    --count;
		}
	    }
	}
	
	if (mfWriteRows(OutputFile,OutputBuffer,1) != 1)
	    return -1;
    }
    return 0;
}



/* ------------------------------------------------------------------
   ReadOrbits() - Read the orbit table file.
   ------------------------------------------------------------------ */

static int ReadOrbits()

{
    FILE *orbfile;

    MESSAGE(1,("Reading orbits from %s",orbname));
    if ((orbfile = sysFopen(orbname,"rb")) == NULL)
	return -1;
    Orbits = imatRead(orbfile);
    if (Orbits == NULL)
    {
	mtxAbort(MTX_HERE,"Error reading orbit table from %s",orbname);
	return -1;
    }
    OrbitSizes = imatRead(orbfile);
    if (OrbitSizes == NULL)
    {
	mtxAbort(MTX_HERE,"Error reading orbit sizes table from %s",orbname);
	return -1;
    }
    Degree = Orbits->Noc;
    NOrbits = OrbitSizes->Noc;
    return 0;
}




/* ------------------------------------------------------------------
   OpenFiles() - Open input and output vector files.
   ------------------------------------------------------------------ */

static int OpenFiles()

{
    /* Open the vector file and allocate row buffer.
       --------------------------------------------- */
    if ((InputFile = mfOpen(vecname)) == NULL)
	return -1;
    if (InputFile->Field < 2)
    {
	mtxAbort(MTX_HERE,"%s: %s",vecname,MTX_ERR_NOTMATRIX);
	return -1;
    }
    if (InputFile->Noc != NOrbits)
    {
	mtxAbort(MTX_HERE,"%s and %s: %s",vecname,orbname,MTX_ERR_INCOMPAT);
	return -1;
    }
    ffSetField(InputFile->Field);
    InputBuffer = ffAlloc(1, NOrbits);

    /* Open the output file and allocate output buffer.
       ------------------------------------------------ */
    MESSAGE(0,("Output is %d x %d\n",InputFile->Nor,Degree));
    OutputFile = mfCreate(resname,ffOrder,InputFile->Nor,Degree);
    if (OutputFile == NULL)
	return -1;
    OutputBuffer = ffAlloc(1, Degree);

    return 0;
}


/* ------------------------------------------------------------------
   Init() - Program initialization.
   ------------------------------------------------------------------ */

static int Init(int argc, char **argv)

{
    App = appAlloc(&AppInfo,argc,argv);
    if (App == NULL)
	return -1;
    if (appGetArguments(App,3,3) < 0)
	return -1;
    vecname = App->ArgV[0];
    orbname = App->ArgV[1];
    resname = App->ArgV[2];

    if (ReadOrbits() != 0)
    {
	mtxAbort(MTX_HERE,"Error reading orbits");
	return -1;
    }
    if (OpenFiles() != 0)
    {
	mtxAbort(MTX_HERE,"Error opening files");
	return -1;
    }
    return 0;
}



static void Cleanup()

{
    if (InputFile != NULL) mfClose(InputFile);
    if (OutputFile != NULL) mfClose(OutputFile);
    if (Orbits != NULL) imatFree(Orbits);
    if (OrbitSizes != NULL) imatFree(OrbitSizes);
    appFree(App);
}


/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char **argv)

{
    int rc;
    if (Init(argc,argv) != 0)
    {
	mtxAbort(MTX_HERE,"Initialization failed");
	return 1;
    }
    rc = uncondense();
    Cleanup();
    return rc;
}




/**
@page prog_zuk zuk - Uncondense Vectors

@section zuk_syntax Command Line
<pre>
zuk @em Options @em Vectors @em Orbits @em Result
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par @em Vectors
  Vectors to uncondense.
@par @em Orbits
  Orbit tables produced by @ref prog_zmo "zmo".
@par @em Result
  Uncondensed vectors.

@section zuk_inp Input Files
@par @em Vectors
  Vectors to uncondense.
@par @em Orbits
  Orbit tables produced by @ref prog_zmo "zmo".

@section zuk_out Output Files
@par @em Result
  Uncondensed vectors.

@section zuk_desc Description
This program reads a matrix which is assumed to be a condensed space of a permutation
representation whose orbits are in the file @em Orbits. The vectors in @em Vectors are
elongated so as to lie in the original permutation space and written out to 
the file @em Result.
@em Orbits must be an orbits file in the format defined by @ref prog_zmo "zmo".
Here is an example:
<pre>
         2 0 4
Space =  1 3 2
         2 0 2

Orbits = (1,2) (3,4,5,6) (7,8,9)

         2 2 0 0 0 0 4 4 4
Result = 1 1 3 3 3 3 2 2 2
         2 2 0 0 0 0 2 2 2
</pre>
*/

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
