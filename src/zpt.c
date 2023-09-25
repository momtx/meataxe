////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Paste matrices.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"
#include <string.h>
#include <stdlib.h>


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static MtxApplicationInfo_t AppInfo = { 
"zpt", "Paste Matrices",
"SYNTAX\n"
"    zpt [-c <NCols>] [-r <NRows>] <Out> [<Inp> ...]\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -c ...................... Set number of columns (matrices only).\n"
"    -r ...................... Set number of rows (matrices only).\n"
"\n"
"ARGUMENTS\n"
"    <Out> ................... Output file name.\n"
"    <Inp> ................... Input file name, '-' to fill with zeroes.\n"
};



static MtxApplication_t *App = NULL;
int nrows = 0;
int ncols = 0;
const char *ofilename;
MtxFile_t* ifile;
MtxFile_t* fileOut;
int fieldOut, norOut, nocOut, maxnor;
int *width, *height;

////////////////////////////////////////////////////////////////////////////////////////////////////

static const char *mkname(int r, int c)
{
    return App->argV[r * ncols + c + 1];
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char **argv)
{   
    App = appAlloc(&AppInfo,argc,argv);

    // Command line options.
    nrows = appGetIntOption(App,"-r",1,1,100);
    ncols = appGetIntOption(App,"-c",1,1,100);

    // Command line arguments.
    if (nrows == 1 && ncols == 1)
    {
	appGetArguments(App,2,1000);
	nrows = App->argC - 1;
    }
    else
    {
    	int names_needed = nrows * ncols + 1;
    	appGetArguments(App,names_needed,names_needed);
    }
    ofilename = App->argV[0];
    width = NALLOC(int,ncols);
    height = NALLOC(int,nrows);

    // Check if we are pasting matrices or permutations.
    MtxFile_t* f = mfOpen(mkname(0,0));
    mfReadHeader(f);
    const uint32_t objectType = mfObjectType(f);
    if (objectType != MTX_TYPE_MATRIX) {
       mtxAbort(MTX_HERE,
             "%s: unsupported object type 0x%lx", mkname(0,0), (unsigned long) objectType);
    }
    mfClose(f);
}




/* ------------------------------------------------------------------
   checksizes()
   ------------------------------------------------------------------ */

static void checksizes()
{
    uint32_t field0 = 0;

    MESSAGE(1,("Checking sizes\n"));
    for (int i = 0; i < nrows; ++i) 
	height[i] = -1;
    for (int k = 0; k < ncols; ++k) 
	width[k] = -1;

    for (int i = 0; i < nrows; ++i)
    {
	for (int k = 0; k < ncols; ++k)
	{
	    const char *c = mkname(i,k);
	    if (!strcmp(c,"-"))
		continue;
            MtxFile_t* f = mfOpen(c);
            mfReadHeader(f);
            const uint32_t objectType = mfObjectType(f);
            if (objectType != MTX_TYPE_MATRIX)
		mtxAbort(MTX_HERE,"%s: %s",c, MTX_ERR_NOTMATRIX);
            const uint32_t field = f->header[0];
            const uint32_t nor = f->header[1];
            const uint32_t noc = f->header[2];
	    mfClose(f);
            if (field0 == 0)
               field0 = field;
	    else if (field0 != field) {
		mtxAbort(MTX_HERE,"%s: wrong field order %lu (expected %lu)",
                      c,(unsigned long) field, (unsigned long) field0);
            }
            if (height[i] < 0)
               height[i] = nor;
            else if (height[i] != nor) {
		mtxAbort(MTX_HERE,"%s: wrong number of rows %lu (expected %lu)",
                      c,(unsigned long) nor, (unsigned long) height[i]);
            }
            if (width[k] < 0)
               width[k] = noc;
            else if (width[k] != noc) {
		mtxAbort(MTX_HERE,"%s: wrong number of columns %lu (expected %lu)",
                      c,(unsigned long) noc, (unsigned long) width[k]);
            }
	}
    }

    fieldOut = field0;
    nocOut = norOut = maxnor = 0;
    for (int i = 0; i < nrows; ++i)
    {
	if (maxnor < height[i])
	    maxnor = height[i];
	norOut += height[i];
    }
    for (int k = 0; k < ncols; ++k)
	nocOut += width[k];

    MESSAGE(0,("Output is %dx%d\n",norOut,nocOut));
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void pasteMatrices()
{   
   ffSetField(fieldOut);
   PTR bufOut = ffAlloc(maxnor, nocOut);
   fileOut = mfCreate(ofilename, fieldOut, norOut, nocOut);
   for (int i = 0; i < nrows; ++i)
   {	
      MESSAGE(1,("Pasting row %d\n",i));
      for (uint32_t row = 0; row < maxnor; ++row)
         ffMulRow(ffGetPtr(bufOut, row, nocOut), FF_ZERO, nocOut);
      int colStart = 0;
      for (int k = 0; k < ncols; ++k)
      {
         const char *c = mkname(i,k);
         if (strcmp(c,"-"))
         {
            MtxFile_t* fileP = mfOpen(c);
            mfReadHeader(fileP);
            const uint32_t norP = fileP->header[1];
            const uint32_t nocP = fileP->header[2];
            PTR piece = ffAlloc(norP, nocP);
            mfReadRows(fileP, piece, norP, nocP);
            mfClose(fileP);
            PTR rowOut = bufOut;
            PTR rowP = piece;
            for (uint32_t r = 0; r < norP; ++r)
            {	
               for (uint32_t c = 0; c < nocP; ++c)
                  ffInsert(rowOut, c+colStart, ffExtract(rowP,c));
               ffStepPtr(&rowP, nocP);
               ffStepPtr(&rowOut, nocOut);
            }
            sysFree(piece);
         }
         colStart += width[k];
      }
      mfWriteRows(fileOut, bufOut, height[i], nocOut);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   init(argc,argv);
   checksizes();
   pasteMatrices();
   appFree(App);
   return 0;
}
 
/**
@page prog_zpt zpt - Paste

@section zpt_syntax Command Line
<pre>
zpt @em Options [-r @em NRows] [-c @em NCols] @em Out @em Inp [@em Inp ...]
</pre>

@par @em Options
Standard options, see @ref prog_stdopts
@par -r @em Rows
Set number of rows.
@par @em Columns
Set number of columns.
@par @em Out
Output file.
@par @em Inp
Input piece.

@section zpt_inp Input Files
@par @em Inp
Input piece.

@section zpt_out Output Files
@par @em Out
Output file.


@section zpt_desc Description
This program reads matrices from one or more input files and pastes the pieces together
to one matrix. The way in which the pieces are put together is controlled by two
parameters, @em NRows and @em NCols. For example,
<pre>
zpt -r 2 -c 3 x aa ab ac ba bb bc
</pre>
would paste together 6 matrices in two rows and three columns.
The resulting matrix is written to "x" and looks like this:
<pre>
aa ab ac
ba bb ba
</pre>
The file name "-" is treated specially: No file is read in,
and the corresponding piece of the output matrix is left empty. This
can be used to calculate the direct sum of two representations. For
example,
<pre>
zpt -r 2 -c 2 A+B A - - B
</pre>
creates the following matrix in block diagonal form:
<pre>
A 0
0 B
</pre>
If only one of @em NRows and @em NCols is specified, the other
parameter is assumed to be one. 
*/

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
