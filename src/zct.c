////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Cut file.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#include <ctype.h>
#include <string.h>
#include <stdlib.h>


#define MAXPIECES 100


/* ------------------------------------------------------------------
   Variables
   ------------------------------------------------------------------ */


static MtxApplicationInfo_t AppInfo = { 
"zct", "Cut Matrices Or Permutations",
"SYNTAX\n"
"    zct <Rows>[:<Columns>] <Input> <Output>\n"
"\n"
"ARGUMENTS\n"
"    <Rows> .................. Rows to cut. A list of integers or ranges\n"
"                              (a-b), separated by commas.\n"
"    <Columns> ............... Columns to cut.\n"
"    <Input> ................. Source file.\n"
"    <Output> ................ Destination file.\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -G ...................... GAP output\n"
"    -m ...................... Calculate the minimal polynomial\n"
"    -f ...................... Factor the polynomial\n"
};

static MtxApplication_t *App = NULL;

static int nRowRanges = 0;
static uint32_t rowlist[MAXPIECES][2];
static int nColRanges = 0;
static uint32_t collist[MAXPIECES][2];

static const char* inputFileName;
static const char* outputFileName;
static uint32_t nor, noc;
static uint32_t onor, onoc;
static MtxFile_t* inputFile;

////////////////////////////////////////////////////////////////////////////////////////////////////

static int list(uint32_t x[MAXPIECES][2], const char **cp)
{   
   const char* c = *cp;
   int n = 0;

   while (isdigit(*c))
   {
      if (n >= MAXPIECES) 
         mtxAbort(MTX_HERE,"Too many row/column ranges (max=%d)", (int)MAXPIECES);
      x[n][0] = (uint32_t) atol(c);
      while (isdigit(*c)) ++c;
      if (*c == '-')
      {	++c;
         if (!isdigit(*c))
         {
            mtxAbort(MTX_HERE,"Invalid range (missing number after '-')");
         }
         x[n][1] = (uint32_t) atol(c);
         while (isdigit(*c)) ++c;
      }
      else
         x[n][1]= x[n][0];
      if (*c == ',') ++c;
      if (x[n][1] < x[n][0]) 
         mtxAbort(MTX_HERE,"Invalid range \"%d-%d\"", x[n][0], x[n][1]);
      ++n;
   }
   *cp = c;
   return n;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int parseargs()

{
    const char *c;

    if (appGetArguments(App,3,3) < 0)
	return -1;

    /* Process the <Range> argument
       ---------------------------- */
    c = App->argV[0];
    nRowRanges = list(rowlist,&c);
    if (*c != 0)
    {
        if (*c != ':' && *c != ';')
	    mtxAbort(MTX_HERE,"Invalid range (':' or ';' expected)");
        ++c;
        nColRanges = list(collist,&c);
    }

    /* Process file names
       ------------------ */
    inputFileName = App->argV[1];
    outputFileName = App->argV[2];

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char **argv)
{
    App = appAlloc(&AppInfo,argc,argv);
    parseargs();

	
    inputFile = mfOpen(inputFileName);
    mfReadHeader(inputFile);
    uint32_t objectType = mfObjectType(inputFile);
	
    if (objectType != MTX_TYPE_MATRIX) {
       mtxAbort(MTX_HERE,"%s: unsupported object type 0x%lx",
             inputFileName, (unsigned long) objectType);
    }
    nor = inputFile->header[1];
    noc = inputFile->header[2];

    if (nColRanges == 0) {	
          nColRanges = 1;
          collist[0][0] = 1;
          collist[0][1] = noc;
    }
    if (nRowRanges == 0) {	
          nRowRanges = 1;
          rowlist[0][0] = 1;
          rowlist[0][1] = nor;
    }
    onor = 0;
    for (uint32_t i = 0; i < nRowRanges; ++i) {
       onor += rowlist[i][1]-rowlist[i][0]+1;
       if (rowlist[i][1] > nor)
          mtxAbort(MTX_HERE, "Row index out of range: %d > %u", rowlist[i][1], (unsigned) nor);
    }
    onoc = 0;
    for (uint32_t i = 0; i < nColRanges; ++i) {
       onoc += collist[i][1]-collist[i][0]+1;
       if (collist[i][1] > nor)
          mtxAbort(MTX_HERE, "Column index out of range: %d > %u", rowlist[i][1], (unsigned)nor);
    }

}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int cutmatrix()
{   
   ffSetField(inputFile->header[0]);
   PTR inputRow = ffAlloc(1, noc);
   PTR outputRow = ffAlloc(1, onoc);
   MtxFile_t* outputFile = mfCreate(outputFileName,inputFile->header[0],onor,onoc);

   for (int i = 0; i < nRowRanges; ++i)
   {	
      const uint32_t row0 = rowlist[i][0]-1;
      sysFseek(inputFile->file, 0);
      mfReadHeader(inputFile);
      sysFseekRelative(inputFile->file, ffRowSizeUsed(noc) * row0);
      for (uint32_t j = rowlist[i][1]-row0; j > 0; --j) {
         mfReadRows(inputFile, inputRow, 1, noc);
         ffMulRow(outputRow, FF_ZERO, onoc);
         uint32_t colOut = 0;
         for (int k = 0; k < nColRanges; ++k) {
            for (uint32_t colIn = collist[k][0]-1; colIn <= collist[k][1]-1; ++colIn) {
               FEL f = ffExtract(inputRow, colIn);
               ffInsert(outputRow, colOut++, f);
            }
         }
         mfWriteRows(outputFile, outputRow, 1, onoc);
      }
   }
   mfClose(outputFile);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   init(argc,argv);
   cutmatrix();
   mfClose(inputFile);
   appFree(App);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// *INDENT-OFF*

/**
@page prog_zct zct - Cut

@section zct_syntax Command Line
<pre>
zct @em Options @em Rows[:@em Columns] @em Input @em Output
</pre>

@par @em Options
Standard options, see @ref prog_stdopts
@par @em Rows
List of rows to cut.
@par @em Columns
List of columns to cut.
@par @em Input
Input matrix.
@par @em Output
Output Matrix.

@section zct_inp Input Files
@par @em Input
Input matrix.

@section zct_out Output Files
@par @em Output
Output matrix.


@section zct_desc Description
This program cuts one or more pieces, specified by @em Rows and @em Columns, out of a matrix
and combines the pieces to a new matrix which is written to the output file.

Both @em Rows and @em Columns are lists of positive integers or ranges (e.g., "13-25")
separated by commas. Ranges can be given in any order, and may be overlapping.
If the columns list is omitted, all columns of the selected rows are cut.

Here are some examples. Assume the input is the following 5 by 10 matrix
<pre>
1 2 3 4 5 6 7 8 9 0
0 1 2 3 4 5 6 7 8 9
0 0 1 2 3 4 5 6 7 8
0 0 0 0 0 0 0 0 0 0
9 8 7 6 5 4 3 2 1 0
</pre>
Then, @b zct would produce the output shown below for different @em Rows:@em Columns lists:
<pre>
Rows:Columns        Result
-----------------   -------------------
1,4-5               1 2 3 4 5 6 7 8 9 0
                    0 0 0 0 0 0 0 0 0 0
                    9 8 7 6 5 4 3 2 1 0

1-3:8-10            8 9 0
                    7 8 9
                    6 7 8

1-2:2,5-7,9         2 5 6 7 9
                    1 4 5 6 8

1-2,5:1-3,9-10      1 2 3 9 0
                    0 1 2 8 9
                    9 8 7 1 0
</pre>
The rows and columns which select the piece need not occur in
ascending order, but the output depends on the ordering.
For example,
<pre>
zct 1,2 output input
zct 2,1 output input
</pre>
both extract the first two rows of `input', but the second form
will also permute the rows. Another example:
<pre>
zct 3-4,1-2:3-4,1-2 inp out
</pre>
would perform the following operation on a 4 by 4 matrix:
<pre>
  inp              out
1 1 2 2          4 4 3 3
1 1 2 2          4 4 3 3
3 3 4 4          2 2 1 1
3 3 4 4          2 2 1 1
</pre>
*/

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
