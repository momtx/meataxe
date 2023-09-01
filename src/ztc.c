////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Trace of a matrix or permutation.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"
#include <stdlib.h>



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */


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

    ffSetField(InputFile->Field); 
    m1 = ffAlloc(1, InputFile->Noc);
    tr = FF_ZERO;
    if ((max = InputFile->Nor) > InputFile->Noc) 
	max = InputFile->Noc;
    for (i = 0; i < max; ++i)
    {
	if (mfReadRows(InputFile,m1,1) != 1)
	{
	    mtxAbort(MTX_HERE,"Cannot read ipnut file");
	    return -1;
	}
	tr = ffAdd(tr,ffExtract(m1,i));
    }
    if (!opt_G)		/* Standard output */
	printf("Trace is %d\n",ffToInt(tr));
    else		/* GAP output */
        printf("MeatAxe.Trace := %s;\n",ffToGap(tr));
    return 0;
}


/* ------------------------------------------------------------------
   trperm() - Trace of a permutation (= number of fixed points)
   ------------------------------------------------------------------ */

static int trperm()
{
    const int degree = InputFile->Nor;

    uint32_t *m1 = NALLOC(uint32_t,degree);
    if (m1 == NULL) 
	return -1;

    mfRead32(InputFile, m1, degree);

    uint32_t tr = 0;
    for (int k = 0; k < degree; ++k)
    {
	if (m1[k] == k) 
	    ++tr;
    }

    if (!opt_G)
	printf("Trace is %ld\n",(unsigned long)tr);
    else
	printf("MeatAxe.Trace := [%ld];\n",(unsigned long)tr);
    return 0;
}


static int Init(int argc, char **argv)

{
    /* Process command line options and arguments.
       ------------------------------------------- */
    if ((App = appAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    opt_G = appGetOption(App,"-G --gap");
    if (opt_G) MtxMessageLevel = -100;
    if (appGetArguments(App,1,1) != 1)
	return -1;
    inpname = App->ArgV[0];

    /* Open the input file.
       -------------------- */
    if ((InputFile = mfOpen(inpname)) == NULL)
    {
	mtxAbort(MTX_HERE,"Error opening input file");
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
	mtxAbort(MTX_HERE,"%s: Unknown type %d",inpname,InputFile->Field);
	return -1;
    }

    return 0;
}


static void Cleanup()

{
    if (App != NULL)
	appFree(App);
    if (InputFile != NULL)
	mfClose(InputFile);
}


/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char **argv)


{
    int rc = 0;

    if (Init(argc,argv) != 0)
    {
	mtxAbort(MTX_HERE,"Initialization failed");
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

@section ztc_syntax Command Line
<pre>
ztc [@em Options] [-G] @em Inp
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par -G
  GAP output.
@par @em Inp
  Input matrix or permutation.

@section ztc_inp Input Files
@par @em Inp
  Input matrix or permutation.

@section ztc_desc Description
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
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
