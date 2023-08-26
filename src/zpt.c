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
"zpt", "Paste Matrices or Permutations",
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
enum { MATRICES, PERMUTATIONS } Mode = MATRICES;
int nrows = 0;
int ncols = 0;
const char *ofilename;
FILE *ifile, *ofile;
int fl, nor, noc, maxnor;
int *width, *height;







/* ------------------------------------------------------------------
   mkname()
   ------------------------------------------------------------------ */

static const char *mkname(int r, int c)

{
    return App->ArgV[r * ncols + c + 1];
}


/* ------------------------------------------------------------------
   Init()
   ------------------------------------------------------------------ */

static int Init(int argc, char **argv)

{   
    MtxFile_t *f;

    App = appAlloc(&AppInfo,argc,argv);
    if (App == NULL)
	return -1;

    /* Command line options.
       --------------------- */
    nrows = appGetIntOption(App,"-r",1,1,100);
    ncols = appGetIntOption(App,"-c",1,1,100);

    /* Command line arguments.
       ----------------------- */ 
    if (nrows == 1 && ncols == 1)
    {
	if (appGetArguments(App,2,1000) < 0)
	    return -1;
	nrows = App->ArgC - 1;
    }
    else
    {
    	int names_needed = nrows * ncols + 1;
    	if (appGetArguments(App,names_needed,names_needed) < 0)
	    return -1;
    }
    ofilename = App->ArgV[0];
    width = NALLOC(int,ncols);
    height = NALLOC(int,nrows);

    /* Check if we are pasting matrices or permutations.
       ------------------------------------------------- */
    if ((f = mfOpen(mkname(0,0))) == NULL)
	return -1;
    if (f->Field < 0)
	Mode = PERMUTATIONS;
    mfClose(f);

    return 0;
}




/* ------------------------------------------------------------------
   checksizes()
   ------------------------------------------------------------------ */

static int checksizes()

{
    int i, k;
    int fl2, nor2, noc2;

    MESSAGE(1,("Checking sizes\n"));
    fl = 0;
    for (i = 0; i < nrows; ++i) 
	height[i] = -1;
    for (k = 0; k < ncols; ++k) 
	width[k] = -1;

    for (i = 0; i < nrows; ++i)
    {
	for (k = 0; k < ncols; ++k)
	{
	    const char *c = mkname(i,k);
	    if (!strcmp(c,"-"))
		continue;
	    if ((ifile = ffReadHeader(c,&fl2,&nor2,&noc2)) == NULL)
		return -1;
	    fclose(ifile);
	    if (fl2 < 2)
	    {
		mtxAbort(MTX_HERE,"%s: %s",c,MTX_ERR_NOTMATRIX);
		return -1;
	    }
	    if (fl == 0)
		fl = fl2;
	    else if (fl != fl2)
	    {
		mtxAbort(MTX_HERE,"%s and %s: %s",mkname(0,0),c,MTX_ERR_INCOMPAT);
		return -1;
	    }
	    if (height[i] == -1) 
		height[i] = nor2;
	    else if (height[i] != nor2)
	    {
		mtxAbort(MTX_HERE,"%s and %s: %s",mkname(i,0),c,MTX_ERR_INCOMPAT);
		return -1;
	    }
	    if (width[k] == -1) 
		width[k] = noc2;
	    else if (width[k] != noc2)
	    {
		mtxAbort(MTX_HERE,"%s and %s: %s",mkname(0,k),c,MTX_ERR_INCOMPAT);
		return -1;
	    }
	}
    }

    /* Calculate nor, noc and maxnor
       ----------------------------- */
    noc = nor = maxnor = 0;
    for (i = 0; i < nrows; ++i)
    {
	if (height[i] == -1) 
	{
	    mtxAbort(MTX_HERE,"Undetermined size");
	    return -1;
	}
	if (height[i] > maxnor)
	    maxnor = height[i];
	nor += height[i];
    }
    for (k = 0; k < ncols; ++k)
    {
	if (width[k] == -1)
	{
	    mtxAbort(MTX_HERE,"Undetermined size");
	    return -1;
	}
	noc += width[k];
    }

    MESSAGE(0,("Output is %dx%d\n",nor,noc));
    return 0;
}



/* ------------------------------------------------------------------
   pastemat() - Paste matrices
   ------------------------------------------------------------------ */

static int pastemat()

{   
    PTR m, piece, x, y;
    int i, k;
    int l, j;
    int pos;

    ffSetField(fl);
    ffSetNoc(noc);
    m = ffAlloc(maxnor, noc);
    if ((ofile = ffWriteHeader(ofilename,fl,nor,noc)) == NULL)
    {
	mtxAbort(MTX_HERE,"Cannot create output file\n");
	return -1;
    }
    for (i = 0; i < nrows; ++i)
    {	
	MESSAGE(1,("Pasting row %d\n",i));
	ffSetNoc(noc);
	x = m;
	for (l = maxnor; l > 0; --l)
	{	
	    ffMulRow(x,FF_ZERO);
	    ffStepPtr(&x, noc);
	}
	pos = 0;
	for (k = 0; k < ncols; ++k)
	{
	    const char *c = mkname(i,k);
	    int fl2, nor2, noc2;

	    if (strcmp(c,"-"))
	    {
		if ((ifile = ffReadHeader(c,&fl2,&nor2,&noc2)) == NULL)
		    return -1;
		ffSetNoc(noc2);
		piece = ffAlloc(nor2, noc2);
		ffReadRows(ifile, piece, nor2, noc2);
		fclose(ifile);
		x = m;
		y = piece;
		ffSetNoc(noc);
 		for (l = 0; l < nor2; ++l)
		{	
		    for (j = 0; j < noc2; ++j)
			ffInsert(x,j+pos,ffExtract(y,j));
		    ffSetNoc(noc2);
		    ffStepPtr(&y, noc2);
		    ffSetNoc(noc);
		    ffStepPtr(&x, noc);
		}
		sysFree(piece);
	    }
	    pos += width[k];
	}
	ffSetNoc(noc);
	ffWriteRows(ofile,m,height[i], noc);
    }
    return 0;
}



static int PastePerms()

{
    MtxFile_t *out, *in;
    int i;
    int degree = 0;
    int nperms = 0;
    long *buf;

    /* Calculate the total number of permutations.
       ------------------------------------------- */
    for (i = 0; i < nrows; ++i)
    {
	if ((in = mfOpen(mkname(i,0))) == NULL)
	    return -1;
	if (in->Field != -1)
	{
	    mtxAbort(MTX_HERE,"%s: %s",mkname(i,0),MTX_ERR_NOTPERM);
	    return -1;
	}
	if (i == 0)
	    degree = in->Nor;
	else if (degree != in->Nor)
	{
	    mtxAbort(MTX_HERE,"Permutations are not compatible");
	    return -1;
	}
	nperms += in->Noc;
	mfClose(in);
    }

    /* Concatenate the permutations.
       ----------------------------- */
    if ((out = mfCreate(ofilename,-1,degree,nperms)) == NULL)
	return -1;
    if ((buf = NALLOC(long,degree)) == NULL)
	return -1;
    for (i = 0; i < nrows; ++i)
    {
	int k;
	if ((in = mfOpen(mkname(i,0))) == NULL)
	    return -1;
	for (k = 0; k < in->Noc; ++k)
	{
	    if (mfReadLong(in,buf,degree) != degree)
    	    {
		mtxAbort(MTX_HERE,"Error reading %s",in->Name);
	    	return -1;
	    }
	    Perm_ConvertOld(buf,degree);
	    if (mfWriteLong(out,buf,degree) != degree)
    	    {
		mtxAbort(MTX_HERE,"Error writing %s",out->Name);
	    	return -1;
	    }
	}
	mfClose(in);
    }
    free(buf);
    mfClose(out);
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
    int rc = 0;
    if (Init(argc,argv) != 0)
    {
	mtxAbort(MTX_HERE,"Initialization failed");
	return 1;
    }
    if (Mode == PERMUTATIONS)
	rc = PastePerms();
    else
    {
	if (checksizes() != 0)
	    return 1;
	if (pastemat() != 0)
	    return 1;
    }
    Cleanup();
    return rc;
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

@subsection ztp_perms Permutations
The program can also paste permutations, i.e., copy permutations from
several files into one file. In this case, "-c" cannot be used. For
example,
<pre>
zpt all perm1 perm2 perm3
</pre>
writes the permutations from "perm1", "perm2", and "perm3" into the
file "all". Each of the input file may contain one or more permutations.
Of course all permutations must have the same degree.
Note: pasting permutations is supported only for compatibility
with older versions of the MeatAxe.

*/

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
