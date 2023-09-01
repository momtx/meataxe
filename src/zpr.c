////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Print a matrix or permutaion.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <ctype.h>
#include <stdlib.h>
#include <string.h>


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */


static long HdrPos;
static uint32_t hdr[3];
static FILE *dest = NULL;   /* Output file */
static const char *inpname = NULL;
static FILE *inpfile = NULL;
static int Gap = 0;		/* -g (GAP mode) */
static int Summary = 0;		/* -s (summary) */



static MtxApplicationInfo_t AppInfo = { 
"zpr", "Print Permutations Or Matrices", 
"SYNTAX\n"
"    zpr [-G] [-s] <Binfile> [<Textfile>]\n"
"\n"
"OPTIONS\n"
"    -G   GAP output\n"
"    -s   Print summary only\n"
"\n"
"FILES\n"
"    <Binfile>   i  A matrix or permutation in binary format\n"
"    <Textfile>  i  The output in text format (default: stdout)\n"
};

static MtxApplication_t *App = NULL;




/* ------------------------------------------------------------------
   PrString(), PrLong() - Prettty printer
   ------------------------------------------------------------------ */

static void PrString(char *c)
{
    static int pos = 0;
    int l = strlen(c);

    if (l + pos >= 78)
    {
	fprintf(dest,"\n");
	pos = 0;
    }
    fprintf(dest,"%s",c);
    for (; *c != 0; ++c)
	if (*c == '\n') pos = 0; else ++pos;
}


static void PrLong(long l)
{
    char tmp[50];
    snprintf(tmp, sizeof(tmp), "%ld",l);
    PrString(tmp);
}


/* ------------------------------------------------------------------
   prmatrix() - Print a matrix in standard format.
   ------------------------------------------------------------------ */

static void prmatrix()
{
    PTR m1;
    FEL f1;
    long loop1, j1;
    int md, mx, iv;

    ffSetField(hdr[0]);
    m1 = ffAlloc(1, hdr[2]);

    if (hdr[0] < 10) {md = 1; mx = 80;}
    else if (hdr[0] < 100) {md = 3; mx = 25;}
    else if (hdr[0] < 1000) {md = 4; mx = 20;}
    else if (hdr[0] < 10000) {md = 5; mx = 15;}
    else {md = 6; mx = 12;}

    fprintf(dest,"matrix field=%ld rows=%ld cols=%ld\n",
          (long)hdr[0],(long)hdr[1],(long)hdr[2]);
    for (loop1 = 1; loop1 <= hdr[1]; ++loop1)
    {	
	ffReadRows(inpfile, m1, 1, hdr[2]);
	iv = 1;
	for (j1 = 0; j1 < hdr[2]; ++j1)
	{
	    f1 = ffExtract(m1,j1);
	    switch (md)
	    {	case 1: fprintf(dest,"%1d",ffToInt(f1)); break;
		case 2: fprintf(dest,"%2d",ffToInt(f1)); break;
		case 3: fprintf(dest,"%3d",ffToInt(f1)); break;
		case 4: fprintf(dest,"%4d",ffToInt(f1)); break;
		case 5: fprintf(dest,"%5d",ffToInt(f1)); break;
		case 6: fprintf(dest,"%6d",ffToInt(f1)); break;
	    }
	    if (iv++ >= mx)
	    {	fprintf(dest,"\n");
		iv = 1;
	    }
	}
	if (iv > 1)
		fprintf(dest,"\n");
    }
}


/* ------------------------------------------------------------------
   prgapmat() - Print a matrix in GAP format.
   ------------------------------------------------------------------ */

static void prgapmat()
{   PTR m1;
   FEL f1;
   FEL gen;
   long loop1, j1;
   int cnt, isprimefield;


   ffSetField(hdr[0]);
   gen = ffGen;		/* Generator */
   isprimefield = (ffChar == ffOrder);

   m1 = ffAlloc(1, hdr[2]);
   PrString("MeatAxe.Matrix := [\n");
   for (loop1 = 1; loop1 <= hdr[1]; ++loop1)
   {	
      ffReadRows(inpfile,m1,1, hdr[2]);
      cnt = 0;
      fprintf(dest,"[");
      for (j1 = 0; j1 < hdr[2]; ++j1)
      {   if (cnt > 75)
         {	fprintf(dest,"\n ");
            cnt = 0;
         }
         f1 = ffExtract(m1,j1);
         if (isprimefield)
         {   FEL f2=FF_ZERO;
            long k=0;
            while (f2 != f1)
            {   f2 = ffAdd(f2,gen);
               ++k;
            }
            fprintf(dest,"%ld",k);
            cnt += k>9999?5:k>999?4:k>99?3:k>9?2:1;
         }
         else
         {   if (f1 == FF_ZERO)
            {   fprintf(dest,"0*Z(%ld)",(long)hdr[0]);
               cnt += 5;
               cnt += hdr[0]>9999?5:hdr[0]>999?4:hdr[0]>99?3:hdr[0]>9?2:1;
            }
            else
            {   FEL f2 = gen;
               long k = 1;
               while (f2 != f1)
               {   f2 = ffMul(f2,gen);
                  ++k;
               }
               fprintf(dest,"Z(%ld)^%ld",(long)hdr[0],k);
               cnt += 4;
               cnt += hdr[0]>9999?5:hdr[0]>999?4:hdr[0]>99?3:hdr[0]>9?2:1;
               cnt += k>9999?5:k>999?4:k>99?3:k>9?2:1;
            }
         }
         if (j1 < hdr[2]-1)
         {	fprintf(dest,",");
            ++cnt;
         }
      }
      fprintf(dest,"]");
      if (loop1 < hdr[1])
         fprintf(dest,",");
      fprintf(dest,"\n");
   }
   fprintf(dest,"]");
   if (isprimefield)
      fprintf(dest,"*Z(%ld)",(long)hdr[0]);
   fprintf(dest,";\n");
}



/* ------------------------------------------------------------------
   prgapimat() - Print an integer matrix in GAP format.
   ------------------------------------------------------------------ */

static void prgapimat()
{
    uint32_t *row;
    long loop1, j1;
    int cnt;

    row = NALLOC(uint32_t,hdr[2]);
    PrString("MeatAxe.Matrix := [\n");
    for (loop1 = 1; loop1 <= hdr[1]; ++loop1)
    {	
	sysRead32(inpfile,row,hdr[2]);
	cnt = 0;
	fprintf(dest,"[");
	for (j1 = 0; j1 < hdr[2]; ++j1)
	{
	    long k = row[j1];
	    if (cnt > 75)
	    {	fprintf(dest,"\n ");
		cnt = 0;
	    }
	    fprintf(dest,"%ld",k);
	    cnt += k>9999?5:k>999?4:k>99?3:k>9?2:1;
	    if (j1 < hdr[2]-1)
	    {	fprintf(dest,",");
		++cnt;
	    }
	}
	fprintf(dest,"]");
	if (loop1 < hdr[1])
	    fprintf(dest,",");
	fprintf(dest,"\n");
    }
    fprintf(dest,"];\n");
}


/* ------------------------------------------------------------------
   prgapperm() - Print a permutation in GAP format.
   ------------------------------------------------------------------ */

static void prgapperm()
{
    long i, pos;
    uint32_t *perm;

    /* Allocate memory for one permutation
      ------------------------------------ */
    perm = NALLOC(uint32_t,hdr[1]);
    if (perm == NULL) 
	mtxAbort(MTX_HERE,"Cannot allocate work space");

    PrString("MeatAxe.Perms := [\n");
    for (pos = 1; pos <= hdr[2]; ++pos)
    {
    	/* Read the next permutation
	   ------------------------- */
	sysRead32(inpfile,perm,hdr[1]);

	/* Print it
	   -------- */
	PrString("    PermList([");
	for (i = 0; i < hdr[1]; ++i)
	{
	    if (i > 0)
	    	PrString(",");
	    PrLong(perm[i] + 1);
	}
	PrString("])");
	if (pos < hdr[2]) PrString(",");
	PrString("\n");
    }
    PrString("];\n");
}



/* ------------------------------------------------------------------
   prgap() - Print a matrix or permutation in GAP format.
   ------------------------------------------------------------------ */

static void prgap()

{
    if (hdr[0] == -1)
	prgapperm();
    else if (hdr[0] == -8)
	prgapimat();
    else if (hdr[0] >= 2)
	prgapmat();
    else
	mtxAbort(MTX_HERE,"Cannot print type %d in GAP format",(int)hdr[0]);
}



/* ------------------------------------------------------------------
   prpol() - Print a polynomial
   ------------------------------------------------------------------ */

static void prpol()

{
    int i;
    Poly_t *p;
    
    sysFseek(inpfile,HdrPos);
    if ((p = polRead(inpfile)) == NULL)
	mtxAbort(MTX_HERE,"Cannot read %s",inpname);
    

    fprintf(dest,"polynomial field=%ld degree=%ld\n",(long)hdr[1],(long)hdr[2]);
    for (i = 0; i <= p->Degree; ++i)
    {
	PrLong(ffToInt(p->Data[i]));
	PrString(" ");
    }
    PrString("\n");
}


/* ------------------------------------------------------------------
   prperm() - Print a permutation in standard format.
   ------------------------------------------------------------------ */

static int prperm()

{
    Perm_t *perm;
    long f1, i;

    sysFseekRelative(inpfile,-3 * 4);
    if ((perm = permRead(inpfile)) == NULL)
    {
	mtxAbort(MTX_HERE,"%s: Cannot read permutation\n");
	return -1;
    }
    
    fprintf(dest,"permutation degree=%d\n",perm->Degree);
    for (i = 0; i < perm->Degree; ++i)
    {
	f1 = perm->Data[i];
	PrLong(f1 + 1);
	PrString(" ");
    }
    PrString("\n");
    return 0;
}

/* ------------------------------------------------------------------
   primat() - Print an integer matrix
   ------------------------------------------------------------------ */

static void primat()

{
    long k, i;
    int32_t *row = NALLOC(int32_t,hdr[2]);
    if (row == NULL) 
	mtxAbort(MTX_HERE,"Cannot allocate work space");

    fprintf(dest,"integer-matrix rows=%ld cols=%ld\n",(long)hdr[1],
	(long)hdr[2]);
    for (i = 0; i < hdr[1]; ++i)
    {
	sysRead32(inpfile,row,hdr[2]);
    	for (k = 0; k < hdr[2]; ++k)
    	{
    	    PrLong(row[k]);
    	    PrString(" ");
    	}
    	PrString("\n");
    }
}



/* ------------------------------------------------------------------
   prmtx() - Print a matrix or permutation in standard format.
   ------------------------------------------------------------------ */

static void prmtx()

{
    if (hdr[0] < 65536U)
	prmatrix();
    else if (hdr[0] == MTX_TYPE_PERMUTATION)
	prperm();
    else if (hdr[0] == MTX_TYPE_POLYNOMIAL)
	prpol();
    else if (hdr[0] == MTX_TYPE_INTMATRIX)
	primat();
    else
	mtxAbort(MTX_HERE,"Cannot print type %d in Mtx format",(int)hdr[0]);
}


static void PrintPermutationSummary()
{
    if (Gap)
    {
	printf("MeatAxe.PermutationCount:=%ld;\n",(long)hdr[2]);
	printf("MeatAxe.PermutationDegree:=%ld;\n",(long)hdr[1]);
    }
    else
    {
	printf("%ld Permutation%s of degree %ld\n",
	    (long)hdr[2], hdr[2] == 1 ? "" : "s", (long)hdr[1]);
    }
}

static void PrintMatrixSummary()
{
    if (Gap)
    {
	printf("MeatAxe.MatrixRows:=%ld;\n",(long)hdr[1]);
	printf("MeatAxe.MatrixCols:=%ld;\n",(long)hdr[2]);
	printf("MeatAxe.MatrixField:=%ld;\n",(long)hdr[0]);
    }
    else
    {
	printf("%ld x %ld matrix over GF(%ld)\n",(long)hdr[1],(long)hdr[2],(long)hdr[0]);
    }
}


static void PrintPolySummary()
{
    if (Gap)
    {
	printf("MeatAxe.PolynomialField=%ld\n",(long)hdr[1]);
	printf("MeatAxe.PolynomialDegree=%ld\n",(long)hdr[2]);
    }
    else
	printf("Polynomial of degree %ld over GF(%ld)\n",(long)hdr[2],(long)hdr[1]);
}


static void PrintImatSummary()
{
    if (Gap)
    {
	printf("MeatAxe.IntegerMatrixRows:=%ld;\n",(long)hdr[1]);
	printf("MeatAxe.IntegerMatrixCols:=%ld;\n",(long)hdr[2]);
    }
    else
	printf("%ld x %ld integer matrix\n",(long)hdr[1],(long)hdr[2]);
}

/* ------------------------------------------------------------------
   PrintSummary() - Print header information
   ------------------------------------------------------------------ */

static void PrintSummary()
{
    if (!Gap)
    	printf("%s: ",inpname);
    if (hdr[0] == -1)
    {
	PrintPermutationSummary();
	sysFseek(inpfile,ftell(inpfile) + hdr[1]*hdr[2]*4);
    }
    else if (hdr[0] >= 2)
    {
	PTR x;
	long l;
	PrintMatrixSummary();
	ffSetField(hdr[0]);
	x = ffAlloc(1, hdr[2]);
	for (l = hdr[1]; l > 0; --l)
	    ffReadRows(inpfile,x,1,hdr[2]);
	free(x);
    }
    else if (hdr[0] == -2)
    {
	PrintPolySummary();
	ffSetField(hdr[1]);
	sysFseek(inpfile,ftell(inpfile) + ffRowSize(hdr[2]+1));
    }
    else if (hdr[0] == -8)
    {
	PrintImatSummary();
	sysFseek(inpfile,ftell(inpfile)+hdr[1]*hdr[2]*4);
    }
    else
	printf("Unknown file format (%ld,%ld,%ld).\n",
              (long)hdr[0],(long)hdr[1],(long)hdr[2]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int readHeader(void)
{
    HdrPos = ftell(inpfile);
    if (!sysTryRead32(inpfile,hdr,3)) return 0;
    if ((hdr[0] > 0x10000 && hdr[0] < MTX_TYPE_INTMATRIX) || hdr[1] < 0 || hdr[2] < 0)
    {
	mtxAbort(MTX_HERE, "%s: %s (%d %d %d)\n",inpname,MTX_ERR_FILEFMT,
	    (int)hdr[0],(int)hdr[1],(int)hdr[2]);
    }
    return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int Init(int argc, char **argv)
{
    if ((App = appAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;

    /* Process command line options
       ---------------------------- */
    Gap = appGetOption(App,"-G --gap");
    Summary = appGetOption(App,"-s --summary");
    if (Gap)
	MtxMessageLevel = -100;	/* Suppress messages in GAP mode */

    /* Process arguments, open files
       ----------------------------- */
    if (appGetArguments(App,1,2) < 0)
	return -1;
    inpname = App->ArgV[0];
    inpfile = sysFopen(inpname,"rb");
    if (inpfile == NULL)
	return -1;
    if (App->ArgC >= 2)
    {  
	dest = sysFopen(App->ArgV[1],"w");
	if (dest == NULL)
	    return -1;
    }
    else
	dest = stdout;
    return 0;
}


/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char **argv)

{
    if (Init(argc,argv) != 0)
    {
	mtxAbort(MTX_HERE,"Initialization failed");
	return 1;
    }

    while (readHeader())
    {
	if (Summary)
	    PrintSummary();
	else if (Gap)
	    prgap();
	else
	    prmtx();
    }
    fclose(inpfile);
    fclose(dest);
    return (EXIT_OK);
}



/**
@page prog_zpr zpr - Print Matrices and Permutations
@see @ref prog_zcv

@section zpr_syntax Command Line
<pre>
zpr [@em Options] [-Gs] @em DataFile [@em TextFile]
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par -G, --gap
  Output in GAP format.
@par -s, --summary
  Show headers only.
@par @em DataFile
  Input file (binary)
@par @em TextFile
  Output file (text)

@section zpr_inp Input Files
@par @em DataFile
  Input file (binary)

@section zpr_out Output Files
@par @em TextFile
  Output file (text)

@section zpr_desc Description
This program prints the contents of a MeatAxe data file in readable
format. The text produced by @b zpr can be converted into binary format by
the @ref prog_zcv "zcv" program.

If there is only one argument on the command line, @b zpr writes to stdout.
A second argument, if present, is taken as the output file name.

To find out the contents of a MeatAxe file, use the -s option. To generate
output readble by GAP, use -G. Both options can be combined.
*/
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
