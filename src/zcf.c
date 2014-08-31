/* ============================= C MeatAxe ==================================
   File:        $Id: zcf.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Change field.
   --------------------------------------------------------------------------
   (C) Copyright 1999 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include <stdlib.h>



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

MTX_DEFINE_FILE_INFO

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
    	MTX_ERROR1("Invalid field order %d",fl2);
	return -1;
    }
    if (fl == -1) 
	return 3;
    if (fl == fl2) 
    {
	MTX_ERROR2("%s is already over GF(%d)",iname,fl);
	return -1;
    }
    else if (fl < fl2)
    {	
	for (f = fl2; f % fl == 0; f /= fl);
	if (f == 1) 
	    return 1;
	else 
	{
	    MTX_ERROR2("Cannot change from GF(%d) to GF(%d)",fl,fl2);
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
	    MTX_ERROR2("Cannot change from GF(%d) to GF(%d)",fl,fl2);
	    return -1;
	}
    }
    return 0;
}




int PermToMat1(const Perm_t *perm, PTR row)

{
    int rc = 0;

    /* Create the output file.
       ----------------------- */
    if ((out = MfCreate(oname,fl2,nor,nor)) == NULL)
	rc = -1;

    /* Convert the permutation.
       ------------------------ */
    if (rc == 0)
    {
	register const long *p = perm->Data;
	int i;
	for (i = 0; rc == 0 && i < nor; ++i)
	{	
	    FfMulRow(row,FF_ZERO);
	    FfInsert(row,p[i],FF_ONE);
	    if (MfWriteRows(out,row,1) != 1)
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
    SysFseek(InputFile,0);
    perm = PermRead(InputFile);
    if (perm == NULL)
	rc = -1;

    /* Allocate workspace.
       ------------------- */
    if (rc == 0)
    {
	FfSetField(fl2);
	FfSetNoc(nor);
	row = FfAlloc(1);
    }

    /* Convert the permutation.
       ------------------------ */
    if (rc == 0)
	rc = PermToMat1(perm, row);

    /* Clean up.
       --------- */
    if (perm != NULL) PermFree(perm);
    if (row != NULL) SysFree(row);

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
	MTX_ERROR("Matrix is too big");
	return -1;
    }
    Buf = NALLOC(FEL,BufNRows * noc);
    if (Buf == NULL)
	return -1;
    FfSetField(fl);
    FfSetNoc(noc);
    if ((RowIn = FfAlloc(1)) == NULL)
	return -1;
    FfSetField(fl2);
    FfSetNoc(noc);
    if ((RowOut = FfAlloc(1)) == NULL)
	return -1;
    return 0;
}



static void FreeBuffer()

{
    if (Buf != NULL)
	SysFree(Buf);
    if (RowIn != NULL)
	SysFree(RowIn);
    if (RowOut != NULL)
	SysFree(RowOut);
}


static int ReadRows(int req)

{
    int to_read;
    int i;
    FEL *tp;

    if ((to_read = req) > BufNRows)
	to_read = BufNRows;
    FfSetField(fl);
    FfSetNoc(noc);
    tp = Buf;
    MESSAGE(1,("Reading %d rows\n",to_read));
    for (i = 0; i < to_read; ++i)
    {
	int k;
	if (FfReadRows(InputFile,RowIn,1) != 1)
	    return -1;
	for (k = 0; k < noc; ++k)
	    *tp++ = FfExtract(RowIn,k);
    }
    return to_read;
}


static int WriteRows(int nrows)

{
    int i;
    FEL *tp;

    FfSetField(fl2);
    FfSetNoc(noc);
    tp = Buf;
    MESSAGE(1,("Writing %d rows\n",nrows));
    for (i = 0; i < nrows; ++i)
    {
	int k;
	for (k = 0; k < noc; ++k)
	    FfInsert(RowOut,k,*tp++);
	if (MfWriteRows(out,RowOut,1) != 1)
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
	if ((out = MfCreate(oname,fl2,nor,noc)) == NULL)
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
	    FfSetField(fl2);
	    for (rp = Buf, k = 0; k < rows_read * noc; ++k, ++rp)
		*rp = FfEmbed(*rp,fl);
	}
	else
	{
	    FfSetField(fl);
	    for (rp = Buf, k = 0; k < rows_read * noc; ++k, ++rp)
		*rp = FfRestrict(*rp,fl2);
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

static int Init(int argc, const char **argv)

{
    /* Process command line options and arguments.
       ------------------------------------------- */
    if ((App = AppAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    if (AppGetArguments(App,3,3) != 3)
	return -1;
    fl2 = atol(App->ArgV[0]);
    iname = App->ArgV[1];
    oname = App->ArgV[2];

    return 0;
}



static void Cleanup()

{
    if (App != NULL)
	AppFree(App);
    if (InputFile != NULL)
	fclose(InputFile);
    if (out != NULL)
	MfClose(out);
}



/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, const char **argv)

{
    int rc = 0;

    if (Init(argc,argv) != 0)
	return -1;

    if ((InputFile = FfReadHeader(iname,&fl,&nor,&noc)) == NULL)
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
	    MTX_ERROR("Conversion failed");
    }

    Cleanup();

    return rc;
}


/**
@page prog_zcf zcf - Change Field

<<<<<<< HEAD
@section syntax Command Line
=======
@section zcf_syntax Command Line
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
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

<<<<<<< HEAD
@section inp Input Files
=======
@section zcf_inp Input Files
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c

@par @em Input
Input file.

<<<<<<< HEAD
@section out Output Files
=======
@section zcf_out Output Files
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c

@par @em Output
Output file.

<<<<<<< HEAD
@section desc Description
=======
@section zcf_desc Description
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c

This program converts between various data types. Currently
there are two kinds of conversions available:
- If @em Input is a matrix, the field is changed to GF(@em q).
  The current field of the input file must be a subfield
  or a superfield of GF(@em q). In the latter case, all matrix
  entries must be in GF(@em q).
- If @em Input is a permutation of degree n, it is converted
  into the corresponding n times n permutation matrix over
  GF(@em q).

<<<<<<< HEAD
@section impl Implementation Details
=======
@section zcf_impl Implementation Details
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
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
