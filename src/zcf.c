////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Change field.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"
#include <stdlib.h>



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */


static const char *iname, *oname;
static FILE *InputFile = NULL;
static MtxFile_t *out = NULL;	    /* Output file */
static int fl;			    /* Field parameter of input file */
static int fl2;			    /* Field to convert to */
static int nor, noc;		    /* Parameters of input file */


static MtxApplicationInfo_t AppInfo = { 
"zcf", "Change Format", 
"SYNTAX\n"
"    zcf " MTX_COMMON_OPTIONS_SYNTAX " <Field> <Input> <Output>\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"\n"
"ARGUMENTS\n"
"    <Field> ................. Desired field order\n"
"    <Input> ................. Input file name\n"
"    <Output> ................ Output file name\n"
};

static MtxApplication_t *App = NULL;





/* ------------------------------------------------------------------
   checkfl() - Check if GF(fl) < GF(fl2)
   ------------------------------------------------------------------ */

static int checkfl()

{	
    int f;
	
    if (fl2 < 2)
    {
    	mtxAbort(MTX_HERE,"Invalid field order %d",fl2);
	return -1;
    }
    if (fl == -1) 
	return 3;
    if (fl == fl2) 
    {
	mtxAbort(MTX_HERE,"%s is already over GF(%d)",iname,fl);
	return -1;
    }
    else if (fl < fl2)
    {	
	for (f = fl2; f % fl == 0; f /= fl);
	if (f == 1) 
	    return 1;
	else 
	{
	    mtxAbort(MTX_HERE,"Cannot change from GF(%d) to GF(%d)",fl,fl2);
	    return -1;
	}
    }
    else if (fl > fl2)
    {	
	for (f = fl; f % fl2 == 0; f /= fl2);
	if (f == 1) 
	    return 1;
	else 
	{
	    mtxAbort(MTX_HERE,"Cannot change from GF(%d) to GF(%d)",fl,fl2);
	    return -1;
	}
    }
    return 0;
}




int permToMat1(const Perm_t *perm, PTR row)

{
    int rc = 0;

    /* Create the output file.
       ----------------------- */
    if ((out = mfCreate(oname,fl2,nor,nor)) == NULL)
	rc = -1;

    /* Convert the permutation.
       ------------------------ */
    if (rc == 0)
    {
	register const long *p = perm->Data;
	int i;
	for (i = 0; rc == 0 && i < nor; ++i)
	{	
	    ffMulRow(row,FF_ZERO);
	    ffInsert(row,p[i],FF_ONE);
	    if (mfWriteRows(out,row,1) != 1)
		rc = -1;
	}
    }
    if (rc == 0)
	MESSAGE(0,("Converted to GF(%d)\n",fl2));
    return 0;
}




/* ------------------------------------------------------------------
   permmat() - Convert permutation -> matrix
   ------------------------------------------------------------------ */

int permmat()

{
    Perm_t *perm = NULL;
    PTR row = NULL;
    int rc = 0;

    /* Read the permutation.
       --------------------- */
    sysFseek(InputFile,0);
    perm = permRead(InputFile);
    if (perm == NULL)
	rc = -1;

    /* Allocate workspace.
       ------------------- */
    if (rc == 0)
    {
	ffSetField(fl2);
	ffSetNoc(nor);
	row = ffAlloc(1, nor);
    }

    /* Convert the permutation.
       ------------------------ */
    if (rc == 0)
	rc = permToMat1(perm, row);

    /* Clean up.
       --------- */
    if (perm != NULL) permFree(perm);
    if (row != NULL) sysFree(row);

    return rc;
}






static int BufNRows = 0;
static FEL *Buf = NULL;
static PTR RowIn = NULL;
static PTR RowOut = NULL;

static int AllocateBuffer()

{
    BufNRows = 1000000 / (sizeof(FEL) * noc);
    if (BufNRows < 1)
    {
	mtxAbort(MTX_HERE,"Matrix is too big");
	return -1;
    }
    Buf = NALLOC(FEL,BufNRows * noc);
    if (Buf == NULL)
	return -1;
    ffSetField(fl);
    ffSetNoc(noc);
    if ((RowIn = ffAlloc(1, noc)) == NULL)
	return -1;
    ffSetField(fl2);
    ffSetNoc(noc);
    if ((RowOut = ffAlloc(1, noc)) == NULL)
	return -1;
    return 0;
}



static void FreeBuffer()

{
    if (Buf != NULL)
	sysFree(Buf);
    if (RowIn != NULL)
	sysFree(RowIn);
    if (RowOut != NULL)
	sysFree(RowOut);
}


static int ReadRows(int req)
{
    int to_read;
    int i;
    FEL *tp;

    if ((to_read = req) > BufNRows)
	to_read = BufNRows;
    ffSetField(fl);
    ffSetNoc(noc);
    tp = Buf;
    MESSAGE(1,("Reading %d rows\n",to_read));
    for (i = 0; i < to_read; ++i)
    {
	int k;
	if (ffReadRows(InputFile,RowIn,1,noc) != 1)
	    return -1;
	for (k = 0; k < noc; ++k)
	    *tp++ = ffExtract(RowIn,k);
    }
    return to_read;
}


static int WriteRows(int nrows)
{
    int i;
    FEL *tp;

    ffSetField(fl2);
    ffSetNoc(noc);
    tp = Buf;
    MESSAGE(1,("Writing %d rows\n",nrows));
    for (i = 0; i < nrows; ++i)
    {
	int k;
	for (k = 0; k < noc; ++k)
	    ffInsert(RowOut,k,*tp++);
	if (mfWriteRows(out,RowOut,1) != 1)
	    return -1;
    }
    return 0;
}




static int ChangeField()
{
    int rc = 0;
    int i;

    /* Allocate work space and create output file.
       ------------------------------------------- */
    rc = AllocateBuffer();
    if (rc == 0)
    {
	if ((out = mfCreate(oname,fl2,nor,noc)) == NULL)
	    rc = -1;
    }

    /* Convert.
       -------- */
    for (i = 0; rc == 0 && i < nor; )
    {
	FEL *rp = Buf;
	int rows_read = ReadRows(nor - i);
	int k;
	if (rows_read <= 0)
	    break;

	/* Convert the marks in <Buf> to GF(<fl2>).
	   ---------------------------------------- */
	MESSAGE(1,("Converting\n"));
	if (fl < fl2)
	{
	    ffSetField(fl2);
	    for (rp = Buf, k = 0; k < rows_read * noc; ++k, ++rp)
		*rp = ffEmbed(*rp,fl);
	}
	else
	{
	    ffSetField(fl);
	    for (rp = Buf, k = 0; k < rows_read * noc; ++k, ++rp)
		*rp = ffRestrict(*rp,fl2);
	}

	if (WriteRows(rows_read) != 0)
	    rc = -1;
	i += rows_read;
    }
    if (i < nor)
	rc = -1;

    FreeBuffer();
    if (rc == 0)
    {
	if (fl < fl2)
	    MESSAGE(0,("Embedded into GF(%d)\n",fl2));
	else
	    MESSAGE(0,("Restricted to GF(%d)\n",fl2));
    }
    return rc;
}

static int Init(int argc, char **argv)

{
    /* Process command line options and arguments.
       ------------------------------------------- */
    if ((App = appAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    if (appGetArguments(App,3,3) != 3)
	return -1;
    fl2 = atol(App->ArgV[0]);
    iname = App->ArgV[1];
    oname = App->ArgV[2];

    return 0;
}



static void Cleanup()

{
    if (App != NULL)
	appFree(App);
    if (InputFile != NULL)
	fclose(InputFile);
    if (out != NULL)
	mfClose(out);
}



/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char **argv)

{
    int rc = 0;

    if (Init(argc,argv) != 0)
	return -1;

    if ((InputFile = ffReadHeader(iname,&fl,&nor,&noc)) == NULL)
	rc = -1;

    /* Convert
       ------- */
    if (rc == 0)
    {
	switch (checkfl())
	{	
	    case 1:	/* Embed into larger field */
		rc = ChangeField();
		break;
	    case 3:
		rc = permmat();
		break;
	}
	if (rc != 0)
	    mtxAbort(MTX_HERE,"Conversion failed");
    }

    Cleanup();

    return rc;
}


/**
@page prog_zcf zcf - Change Field

@section zcf_syntax Command Line
<pre>
zcf @em Options @em q @em Input @em Output
</pre>

@par @em Options
Standard options, see @ref prog_stdopts

@par @em q
New field order.

@par @em Input
Input file.

@par @em Output
Output file.

@section zcf_inp Input Files

@par @em Input
Input file.

@section zcf_out Output Files

@par @em Output
Output file.

@section zcf_desc Description

This program converts between various data types. Currently
there are two kinds of conversions available:
- If @em Input is a matrix, the field is changed to GF(@em q).
  The current field of the input file must be a subfield
  or a superfield of GF(@em q). In the latter case, all matrix
  entries must be in GF(@em q).
- If @em Input is a permutation of degree n, it is converted
  into the corresponding n times n permutation matrix over
  GF(@em q).

@section zcf_impl Implementation Details
For matrices, the conversion is done in two steps. First,
all entries of the matrix are converted to integers.
Then, they are mapped to the new field and reassembled into rows.
The result is written out row by row.

In case of permutations the output matrix is generated row by
row by inserting ones at the positions specified by the
permutation.

If the input is a matrix, the whole matrix must fit into memory.
Additionally, the program needs n⋅m⋅s bytes of memory, where m and n
are the dimensions of the input matrix and s=1 for the small
arithmetic version and s=2 for the big version.
In case of permutations, the input permutation and one row
of the output file must fit into memory.
*/
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
