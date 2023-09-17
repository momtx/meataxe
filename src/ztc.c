////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Trace of a matrix or permutation.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"
#include <stdlib.h>


static int opt_G = 0;		/* GAP output */
static MtxFile_t *InputFile = NULL;
static const char *inpname = NULL;
static MtxApplication_t *App = NULL;

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

////////////////////////////////////////////////////////////////////////////////////////////////////

static void trmat()
{
    ffSetField(InputFile->header[0]); 
    const uint32_t nor = InputFile->header[1];
    const uint32_t noc = InputFile->header[2];
    const uint32_t min = nor < noc ? nor : noc;
    PTR m1 = ffAlloc(1, noc);
    FEL tr = FF_ZERO;
    for (uint32_t i = 0; i < min; ++i)
    {
	mfReadRows(InputFile,m1,1,noc);
	tr = ffAdd(tr,ffExtract(m1,i));
    }
    if (!opt_G)
	printf("Trace is %lu\n",(unsigned long) ffToInt(tr));
    else
        printf("MeatAxe.Trace := %s;\n",ffToGap(tr));
    sysFree(m1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int trperm()
{
    const uint32_t degree = InputFile->header[1];
    uint32_t *m1 = NALLOC(uint32_t,degree);
    mfRead32(InputFile, m1, degree);
    uint32_t tr = 0;
    for (int k = 0; k < degree; ++k)
    {
	if (m1[k] == k) 
	    ++tr;
    }
    sysFree(m1);

    if (!opt_G)
	printf("Trace is %lu\n",(unsigned long)tr);
    else
	printf("MeatAxe.Trace := [%lu];\n",(unsigned long)tr);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char **argv)
{
    App = appAlloc(&AppInfo,argc,argv);
    opt_G = appGetOption(App,"-G --gap");
    if (opt_G) MtxMessageLevel = -100;
    appGetArguments(App,1,1);
    inpname = App->ArgV[0];
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   init(argc,argv);
   InputFile = mfOpen(inpname);
   mfReadHeader(InputFile);
   const uint32_t objectType = mfObjectType(InputFile);
   if (objectType == MTX_TYPE_MATRIX)
      trmat();
   else if (objectType == MTX_TYPE_PERMUTATION)
      trperm();
   else
      mtxAbort(MTX_HERE,"%s: Unknown object type 0x%lx",inpname,(unsigned long)objectType);
   appFree(App);
   mfClose(InputFile);
   return 0;
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
