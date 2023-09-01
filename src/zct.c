////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Cut file.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#include <ctype.h>
#include <string.h>
#include <stdlib.h>


#define MAXPIECES 10


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

static int nrows = 0;
static long rowlist[MAXPIECES][2];
static int ncols = 0;
static long collist[MAXPIECES][2];
static const char *savedpos;
static const char *ifilename;
static const char *ofilename;
static int fl, nor, noc;
static int onor, onoc;
static FILE *InputFile, *OutputFile;



/* ------------------------------------------------------------------
   err() - Print error message and exit
   ------------------------------------------------------------------ */

static void err(int c)

{   
    fprintf(stderr,"ZCT ERROR - ");
    switch (c)
    {	
	case 'a':
	    fprintf(stderr,"BAD RANGE\n");
	    break;
	case 'x':
	    fprintf(stderr,"TOO MANY PIECES\n");
	    break;
	case 'P':
	    fprintf(stderr,"CANNOT CUT COLUMNS FROM A PERMUTATION\n");
	    break;
	case 'o':
	    fprintf(stderr,"CANNOT CUT THE REQUESTED PIECE\n");
	    break;
    }
    exit(EXIT_ERR);
}


/* ------------------------------------------------------------------
   list() - Scan a list of rows/columns
   ------------------------------------------------------------------ */

static int list(long x[MAXPIECES][2], const char *c)

{   
    int n = 0;

    while (isdigit(*c))
    {
	if (n >= MAXPIECES) err('x');
	x[n][0] = atol(c);
	while (isdigit(*c)) ++c;
	if (*c == '-')
	{	++c;
		if (!isdigit(*c))
		{
    	    	    mtxAbort(MTX_HERE,"Invalid range (missing number after '-')");
		    return -1;
		}
		x[n][1] = atol(c);
		while (isdigit(*c)) ++c;
	}
	else
		x[n][1]= x[n][0];
	if (*c == ',') ++c;
	if (x[n][1] < x[n][0]) err('a');
	++n;
    }
    savedpos = c;
    return n;
}




/* ------------------------------------------------------------------
   parseargs()
   ------------------------------------------------------------------ */

static int parseargs()

{
    const char *c;

    if (appGetArguments(App,3,3) < 0)
	return -1;

    /* Process the <Range> argument
       ---------------------------- */
    c = App->ArgV[0];
    nrows = list(rowlist,c);
    c = savedpos;
    if (*c != 0)
    {
        if (*c != ':' && *c != ';')
	{
	    mtxAbort(MTX_HERE,"Invalid range (':' or ';' expected)");
	    return -1;
	}
        ++c;
        ncols = list(collist,c);
    }

    /* Process file names
       ------------------ */
    ifilename = App->ArgV[1];
    ofilename = App->ArgV[2];

    return 0;
}



/* ------------------------------------------------------------------
   init()
   ------------------------------------------------------------------ */

static int init()

{
    int i;
	
    if ((InputFile = ffReadHeader(ifilename,&fl,&nor,&noc)) == NULL)
	return -1;
    if (fl == -1)	/* Is it a permutation? */
    {	
	if (ncols > 0) err('P');
    }
    else
    {	
	if (ncols == 0)
	{	
	    ncols = 1;
	    collist[0][0] = 1;
	    collist[0][1] = noc;
	}
    }

    if (nrows == 0)
	onor = nor;
    else
    {	
	onor = 0;
	for (i = 0; i < nrows; ++i)
	    onor += rowlist[i][1]-rowlist[i][0]+1;
    }
    if (ncols == 0)
	onoc = noc;
    else
    {	
	onoc = 0;
	for (i = 0; i < ncols; ++i)
	    onoc += collist[i][1]-collist[i][0]+1;
    }

    return 0;
}




static int Init(int argc, char **argv)

{
    if ((App = appAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    if (parseargs() != 0)
	return -1;
    if (init() != 0)
	return -1;
    return 0;
}






/* ------------------------------------------------------------------
   cutperm()
   ------------------------------------------------------------------ */

static void cutperm()
{
    int i;
    uint32_t *x, *y;

    if (rowlist[nrows-1][1] > noc) err('o');

    // Read permutations
    y = x = NALLOC(uint32_t, nor * onor);
    for (i = 0; i < nrows; ++i)
    {
	int n = nor * (rowlist[i][1]-rowlist[i][0]+1);
	sysFseek(InputFile,12 + (rowlist[i][0] - 1) * nor * 4);
	sysRead32(InputFile,y,n);
	y += n;
    }

    // Write output
    OutputFile = ffWriteHeader(ofilename,(long)-1,nor,onor);
    sysWrite32(OutputFile, x, onor*nor);
}


/* ------------------------------------------------------------------
   cutmatrix()
   ------------------------------------------------------------------ */

static int cutmatrix()

{   int i, ii;
    int k, kk, pos;
    PTR x,y,row;

    if (rowlist[nrows-1][1] > nor) err('o');
    if (collist[ncols-1][1] > noc) err('o');
    ffSetField(fl);
    row = ffAlloc(1, noc);
    x = ffAlloc(onor, onoc);
    y = x;
    for (i = 0; i < nrows; ++i)
    {	
	ffSeekRow(InputFile, rowlist[i][0]-1, noc);
	for (k = 0; k <= rowlist[i][1]-rowlist[i][0]; ++k)
	{   
	    ffReadRows(InputFile,row,1, noc);
	    pos = 0;
	    ffMulRow(y,FF_ZERO, onoc);
	    for (ii = 0; ii < ncols; ++ii)
	    {
		for (kk = collist[ii][0]-1; kk < collist[ii][1]; ++kk)
		{   FEL f = ffExtract(row,kk);
		    ffInsert(y,pos,f);
		    ++pos;
		}
	    }
	    ffStepPtr(&y, onoc);
	}
    }

    /* Write output
       ------------ */
    if ((OutputFile = ffWriteHeader(ofilename,fl,onor,onoc)) == NULL)
	return -1;
    ffWriteRows(OutputFile,x,onor, onoc);
    return 0;
}

static void Cleanup()

{
    if (OutputFile != NULL)
	fclose(OutputFile);
    if (InputFile != NULL)
	fclose(InputFile);
    appFree(App);
}


int main(int argc, char **argv)

{
    int rc = 0;

    if (Init(argc,argv) != 0)
    {
	mtxAbort(MTX_HERE,"Initialization failed");
	return -1;
    }
    if (fl == -1)
	cutperm();
    else
	rc = cutmatrix();

    Cleanup();
    return rc;
}



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
Input matrix or permutation.
@par @em Output
Output file.

@section zct_inp Input Files
@par @em Input
Input matrix or permutations.

@section zct_out Output Files
@par @em Output
Output file.


@section zct_desc Description
This program cuts a piece, specified by @em Rows and @em Columns,
out of the file @em Input, and writes the piece to @em Output.
The input may be a matrix or a set of permutations.
Both @em Rows and @em Columns are lists of positive
integers or ranges (e.g., "13-25") separated by commas.
If the input is a matrix, the corresponding rows and columns are cut,
and the resulting rectangular pieces
are combined into one rectangular matrix which is written to the
output file.
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
With permutations the program works in the same way as with
matrices.  Each permutation is treated as a row. The @em Columns
list must be empty in  this case, because @b zct can cut only entire
permutations.


@section zct_impl Implementation Details
The number of entries in the @em Rows and @em Columns list
must not be greater than 10. One row (or permutation, respectively)
of the input file and the whole result of the cut must fit into
memory.


*/
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
