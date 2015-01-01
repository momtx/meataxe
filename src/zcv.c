////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Convert a matrix or permutation from ASCII (readable)
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////


#define MAXLINE 4000	/* Max. input line size */

#include "meataxe.h"
#include <ctype.h>
#include <string.h>
#include <stdlib.h>




/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

MTX_DEFINE_FILE_INFO


static int GrpLibFormat = 0;	/* File is in Group Library Format */
static int fl;
static int mod;
static int nor,noc;
static FILE *src = NULL;		/* Input file */
static FILE *out;			/* Output file */
static char lbuf[MAXLINE] = {0};	/* Input line buffer */
static char *lptr = lbuf;		/* Read pointer */
static const char *inpname = "[stdin]";
static const char *outname = "";
static int MemberCount = 0;		/* Number of members */

static MtxApplicationInfo_t AppInfo = { 
"zcv","Convert Text to Binary Format", 
"SYNTAX\n"
"    zcv <Inp> <Out>\n"
"\n"
"ARGUMENTS\n"
"    <Inp> ................... Input file (MeatAxe text format). '-' for stdin\n"
"    <Out> ................... Output file (MeatAxe binary format)\n"
"\n"
"FILES\n"
"    <Inp> ................... I Text file\n"
"    <Out> ................... O Binary file\n"
};


static MtxApplication_t *App = NULL;






/* ------------------------------------------------------------------
   readline() - Read next input line, skip empty lines and strip
   	comments. Returns 0 on success, 1 on end-of-file.
   ------------------------------------------------------------------ */

static int readline()

{
    char *c;
    int mt = 1;

    while (mt)
    {	lbuf[0] = 0;
	if (feof(src)) return 1;
	fgets(lbuf,sizeof(lbuf),src);
	if (ferror(src)) 
	{
	    MTX_ERROR("Unexpected end of input file");
	    return -1;
	}
	for (c = lbuf; *c != 0 && *c != '#'; ++c)
	    if (!isspace(*c)) mt = 0;
	*c = 0;
    }
    lptr = lbuf;
    return 0;
}


/* ------------------------------------------------------------------
   readlong() - Read integer number
   ------------------------------------------------------------------ */

static long readlong()

{
    long l;
    int minus = 0;

    while (!isdigit(*lptr) && *lptr != '-')
    {
	while (*lptr != 0 && !isdigit(*lptr) && *lptr != '-') ++lptr;
	if (*lptr == 0)
	{
	    if (readline()) 
	    {
		MTX_ERROR("Unexpected end of input file");
		return -1;
	    }
	}
    }
    if (*lptr == '-') { minus = 1; ++lptr; }
    if (!isdigit(*lptr))
    {
	MTX_ERROR1("%s: Bad file format",inpname);
	return -1;
    }
    for (l = 0; isdigit(*lptr); )
    {	
	l *= 10;
	switch (*lptr)
	{   case '0': break;
	    case '1': l += 1; break;
	    case '2': l += 2; break;
	    case '3': l += 3; break;
	    case '4': l += 4; break;
	    case '5': l += 5; break;
	    case '6': l += 6; break;
	    case '7': l += 7; break;
	    case '8': l += 8; break;
	    case '9': l += 9; break;
	}
	++lptr;
	if (mod != 2 && mod != 5 && fl >= 2 && fl <= 9) break;
    }
    return minus ? -l : l;
}


/* ------------------------------------------------------------------
   WriteHeader() - Write a header consisting of three integers
   ------------------------------------------------------------------ */

static void WriteHeader(long a, long b, long c)

{
    long hdr[3];
    hdr[0] = a;
    hdr[1] = b;
    hdr[2] = c;
    if (SysWriteLong(out,hdr,3) != 3)
	MTX_ERROR("Cannot write header");
}


/* ------------------------------------------------------------------
   convmatrix() - Convert matrix in fixed format (mode 1)
   ------------------------------------------------------------------ */

static void convmatrix()

{	
    int i, j;
    PTR m1;
    long val = 0;
    int inp;

    if (fl > 9) 
    {
	MTX_ERROR1("Mode 1 not allowed for GF(%d)",fl);
	return;
    }
    FfSetField(fl); 
    FfSetNoc(noc);
    m1 = FfAlloc(1);
    WriteHeader(fl,nor,noc);
    MESSAGE(0,("%dx%d matrix over GF(%d)\n",nor,noc,fl));
    for (i = 1; i <= nor; ++i)
    {	FfMulRow(m1,FF_ZERO);
	inp = 81;
	for (j = 0; j < noc; ++j)
	{	if (inp >= 80)	/* read next line */
		{	
		    memset(lbuf,0,sizeof(lbuf));
		    if (readline()) 
			return;
		    inp = 0;
		}
		switch (lbuf[inp++])
		{	case ' ':
			case '0': val = 0; break;
			case '1': val = 1; break;
			case '2': val = 2; break;
			case '3': val = 3; break;
			case '4': val = 4; break;
			case '5': val = 5; break;
			case '6': val = 6; break;
			case '7': val = 7; break;
			case '8': val = 8; break;
			default: 
			    MTX_ERROR1("%s: Bad file format (Digit expected)",inpname);
		}
		if (val > fl) 
		    MTX_ERROR1("%s: Bad file format",inpname);
		FfInsert(m1,j,FfFromInt(val));
	}
	FfWriteRows(out,m1,1);
    }
}


/* ------------------------------------------------------------------
   conv23456() - Convert matrix in free format (mode 3,4,5,6)
   ------------------------------------------------------------------ */

static void conv3456()

{	
    int i, j;
    PTR m1;
    long val;

    MESSAGE(0,("%dx%d matrix over GF(%d)\n",nor,noc,fl));
    FfSetField(fl); 
    FfSetNoc(noc);
    m1 = FfAlloc((long)1);
    WriteHeader(fl,nor,noc);
    for (i = 1; i <= nor; ++i)
    {	FfMulRow(m1,FF_ZERO);
	for (j = 0; j < noc; ++j)
	{
	    val = readlong();
	    if (mod == 5)
	    {	val %= FfChar;
		if (val < 0) val += FfChar;
	    }
	    FfInsert(m1,j,FfFromInt(val));
	}
	FfWriteRows(out,m1,1);
    }
}


/* ------------------------------------------------------------------
   ConvertMatrix() - Convert matrix (new format)
   ------------------------------------------------------------------ */

static void ConvertMatrix()

{
    long i, j;
    PTR m1;
    long val;

    MESSAGE(0,("%dx%d matrix over GF(%d)\n",nor,noc,fl));
    FfSetField(fl); 
    FfSetNoc(noc);
    m1 = FfAlloc((long)1);
    WriteHeader(fl,nor,noc);
    for (i = 1; i <= nor; ++i)
    {
    	FfMulRow(m1,FF_ZERO);
	for (j = 0; j < noc; ++j)
	{
	    val = readlong();
	    FfInsert(m1,j,FfFromInt(val));
	}
	FfWriteRows(out,m1,1);
    }
}


/* ------------------------------------------------------------------
   ConvertIntMatrix() - Convert integer matrix
   ------------------------------------------------------------------ */

static void ConvertIntMatrix()

{
    long i, j;
    long *x;

    MESSAGE(0,("%dx%d integer matrix\n",nor,noc));
    x = NALLOC(long,noc);
    WriteHeader(-8,nor,noc);	    /* 8 = T_IMAT */
    for (i = 1; i <= nor; ++i)
    {
	for (j = 0; j < noc; ++j)
	    x[j] = readlong();
	SysWriteLong(out,x,noc);
    }
    free(x);
}


/* ------------------------------------------------------------------
   ConvertPermutation() - Convert permutation (new format)
   ------------------------------------------------------------------ */

static void ConvertPermutation()

{
    long i;
    long *buf;
    long kk;


    MESSAGE(0,("Permutation on %d points\n",nor));
    buf = NALLOC(long,nor);
    if (buf == NULL) 
	MTX_ERROR("Cannot allocate permutation: %S");
    WriteHeader(-1,nor,1);

    for (i = 0; i < nor; ++i)
    {
	kk = readlong();
	buf[i] = kk - 1;
	if (kk < 1 || kk > nor)
	{
	    MTX_ERROR3("%s: Invalid point %d in permutation of degree %d",
		inpname,(int)kk,nor);
	}
    }
    if (SysWriteLong(out,buf,nor) != nor)
	MTX_ERROR1("Cannot write to %s",outname);
}


/* ------------------------------------------------------------------
   ConvertPolynomial() - Convert polynomial (new format)
   ------------------------------------------------------------------ */

static void ConvertPolynomial()

{
    long i;
    Poly_t *p;

    MESSAGE(0,("Polynomial of degree %d over GF(%d)\n",nor,fl));
    FfSetNoc(fl);
    p = PolAlloc(fl,nor);
    if ((out = SysFopen(outname,FM_CREATE)) == NULL)
	MTX_ERROR1("Cannot open %s: %S",outname);
    for (i = 0; i <= nor; ++i)
    {
	long kk = readlong();
	p->Data[i] = FfFromInt(kk);
    }
    PolWrite(p,out);
}


/* ------------------------------------------------------------------
   convperm() - Convert permutation (mode 2)
   ------------------------------------------------------------------ */

void convperm()		/* mode 2 */

{
    long i, val;
    PTR m1;

    MESSAGE(0,("%dx%d permutation matrix over GF(%d)\n",nor,noc,fl));
    FfSetField(fl); 
    FfSetNoc(noc); 
    m1 = FfAlloc((long)1);
    WriteHeader(fl,nor,noc);
    for (i = 1; i <= nor; ++i)
    {	
	val = readlong();
	FfMulRow(m1,FF_ZERO);
	FfInsert(m1,val - 1,FF_ONE);
	FfWriteRows(out,m1,1);
    }
}


/* ------------------------------------------------------------------
   conv1213() - Convert permutation (mode 12, 13)
   ------------------------------------------------------------------ */

void conv1213()		/* modes 12, 13 */

{
    long nper, i;
    long *buf;
    long kk = 0, f1, f2;


    MESSAGE(0,("Permutation on %d points\n",nor));
    buf = NALLOC(long,nor);
    if (buf == NULL) 
    {
	MTX_ERROR("Cannot allocate permutation: %S");
	return;
    }
    WriteHeader(-fl,nor,noc);
    for (nper = noc; nper != 0; --nper)
    {	
	for (i = 0; i < nor; ++i)
	{
	    switch (mod)
	    {	case 12:
			kk = readlong();
			break;
		case 13:
			f1 = readlong();
			f2 = readlong();
			kk = (f1-1)*fl + f2 + 1;
			break;
	    }
	    buf[i] = kk;
	}
	if (SysWriteLong(out,buf,nor) != nor)
	    MTX_ERROR1("Cannot write %s: %S",outname);
    }
}


/* ------------------------------------------------------------------
   ReadHeader() - Read the next header. Returns 0 on end of file,
   1 on success.
   ------------------------------------------------------------------ */

static int ReadHeader(void)

{
    if (readline())
    	return 0;
    return 1;
}


/* ------------------------------------------------------------------
   Convert() - Convert one member.
   ------------------------------------------------------------------ */

static void Convert(void)

{
    static char sfl[20], smode[20], snor[20], snoc[20];
    char *c;

    /* Check for new file format
       ------------------------- */

    if ((c = strstr(lbuf,"MeatAxeFileInfo")) != NULL)
    {
	char *d = lbuf;
	for (c += 15; *c != 0 && *c != '"'; ++c);
	if (*c != '"') 
	    MTX_ERROR1("%s: Bad file format",inpname);
	for (++c; *c != '"' && *c != 0; ++c)
	    *d++ = *c;
	*d = 0;
	GrpLibFormat = 1;
    }

    if (!strncmp(lbuf,"matrix",6))
    {
    	char *c;
    	fl = nor = noc = -1;
    	for (c=strtok(lbuf+6," \t\n"); c!=NULL; c=strtok(NULL," \t\n"))
    	{
    	    if (!strncmp(c,"field=",6)) fl = atol(c+6);
    	    else if (!strncmp(c,"nor=",4)) nor = atol(c+4);
    	    else if (!strncmp(c,"rows=",5)) nor = atol(c+5);
    	    else if (!strncmp(c,"noc=",4)) noc = atol(c+4);
    	    else if (!strncmp(c,"cols=",5)) noc = atol(c+5);
    	    else MTX_ERROR1("%s: Bad file format",inpname);
    	}
    	if (nor == -1 || noc == -1 || fl == -1) MTX_ERROR1("%s: Bad file format",inpname);
    	readline();
    	ConvertMatrix();
    	return;
    }
    if (   !strncmp(lbuf,"integer matrix",14)
        || !strncmp(lbuf,"integer-matrix",14))
    {
    	char *c;
    	fl = nor = noc = -1;
    	for (c=strtok(lbuf+14," \t\n"); c!=NULL; c=strtok(NULL," \t\n"))
    	{
    	    if (!strncmp(c,"nor=",4)) nor = atol(c+4);
    	    else if (!strncmp(c,"rows=",5)) nor = atol(c+5);
    	    else if (!strncmp(c,"noc=",4)) noc = atol(c+4);
    	    else if (!strncmp(c,"cols=",5)) noc = atol(c+5);
    	    else MTX_ERROR1("%s: Bad file format",inpname);
    	}
	if (nor == -1 || noc == -1) MTX_ERROR1("%s: Bad header format",inpname);
    	readline();
    	ConvertIntMatrix();
    	return;
    }
    else if (!strncmp(lbuf,"permutation",11))
    {
    	char *c;
    	fl = nor = -1;
    	noc = 1;
    	for (c=strtok(lbuf+11," \t\n"); c!=NULL; c=strtok(NULL," \t\n"))
    	{
    	    if (!strncmp(c,"degree=",7)) nor = atol(c+7);
    	    else if (!strncmp(c,"deg=",4)) nor = atol(c+4);
    	    else MTX_ERROR1("%s: Bad header format",inpname);
    	}
    	if (nor == -1) MTX_ERROR1("%s: Bad header format",inpname);
    	readline();
    	ConvertPermutation();
    	return;
    }
    else if (!strncmp(lbuf,"polynomial",10))
    {
    	char *c;
    	fl = nor = -1;
    	noc = 1;
    	for (c=strtok(lbuf+11," \t\n"); c!=NULL; c=strtok(NULL," \t\n"))
    	{
    	    if (!strncmp(c,"degree=",7)) nor = atol(c+7);
    	    else if (!strncmp(c,"deg=",4)) nor = atol(c+4);
    	    else if (!strncmp(c,"field=",6)) fl = atol(c+6);
    	    else MTX_ERROR1("%s: Bad polynomial header format",inpname);
    	}
    	if (nor < 0 || fl < 2) 
	    MTX_ERROR3("%s: Bad header: fl=%d, deg=%d",inpname,fl,nor);
    	readline();
    	ConvertPolynomial();
    	return;
    }

    /* Use old file format
       ------------------- */
    strncpy(smode,lbuf,2); smode[2] = 0;
    strncpy(sfl,lbuf+2,6); sfl[6] = 0;
    strncpy(snor,lbuf+8,6); snor[6] = 0;
    strncpy(snoc,lbuf+14,6); snoc[6] = 0;
    if (sscanf(lbuf,"%d%d%d%d",&mod,&fl,&nor,&noc) != 4)
    {
        sscanf(smode,"%d",&mod);
        sscanf(sfl,"%d",&fl);
        sscanf(snor,"%d",&nor);
        sscanf(snoc,"%d",&noc);
    }

    if (mod != 1) readline();
    switch (mod)
    {
	case 1:
	    convmatrix();
	    break;
	case 2:
	    convperm();
	    break;
	case 3:
	case 4:
	case 5:
	case 6:
	    conv3456();
	    break;
	case 12:
	case 13:
	    conv1213();
	    break;
	default:
	    MTX_ERROR2("%s: Unknown mode %d",inpname,mod);
    }
}





static int Init(int argc, const char **argv)

{	
    if ((App = AppAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    if (AppGetArguments(App,2,2) < 0)
	return -1;

    /* Open input file.
       ---------------- */
    inpname = App->ArgV[0];
    if (strcmp(inpname,"-"))
    {
	src = SysFopen(inpname,FM_READ|FM_TEXT|FM_LIB);
	if (src == NULL)
	{
	    MTX_ERROR1("Cannot open %s",inpname);
	    return -1;
	}
    }
    else
	src = stdin;

    /* Open output file.
       ----------------- */
    outname = App->ArgV[1];
    out = SysFopen(outname,FM_CREATE);
    if (out == NULL)
    {
	MTX_ERROR1("Cannot open %s for output",outname);
	return -1;
    }

    return 0;
}


static void Cleanup()

{
    fclose(out);
    if (App != NULL) AppFree(App);
}




/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, const char **argv)

{	
    if (Init(argc,argv) != 0)
	return 0;

    while (ReadHeader())
    {
    	Convert();
	++MemberCount;
    }
    if (MemberCount == 0)
    	MESSAGE(0,("Warning: %s is empty",inpname));

    Cleanup();
    return 0;
}




/**
@page prog_zcv zcv - Convert Text to Binary Format
@see @ref prog_zpr

@section zcv_syntax Command Line
<pre>
zcv @em Options @em TextFile @em DataFile
</pre>

@par @em Options
Standard options, see @ref prog_stdopts

@par @em TextFile
Input file (text)

@par @em DataFile
Output file (binary)

@section zcv_inp Input Files

@par @em TextFile
Input file (text)

@section zcv_out Output Files

@par @em DataFile
Output file (binary)

@section zcv_desc Description
This program converts a text file into binary format.
If the input file name is "-", input is read from stdin. 

@par Text File Format
The text file is interpreted line by line. Empty lines are ignored 
completetly. If the line contains one or more "#", which may occur at 
any position, the first "#" and all remaining characters on this line
are ingored.

Each object, for example a matrix or a permutation, consists of a one 
line header followed by the data. Here are some examples of possible
header formats:
<pre>
matrix field=16 rows=10 cols=10
permutation degree=10026
polynomial field=2 degree=23
</pre>
The header may be given in a different format, for example
<pre>
MeatAxeFileInfo := "matrix field=5 rows=100 cols=100";
</pre>
Other header formats are supported for compatibility with older versions 
of the MeatAxe.

After an old type 1 header, the format of the data follwing the header is 
fixed. There must be at most 80 characters per line, lines must be filled 
as much as possible, and each row of the matrix starts with a new input 
line. There are no blanks allowed to separate numbers. In all other modes, 
the data part consists of a sequence of integers in free format, separated 
by any combination of blanks, tabs, or newlines.


*/
