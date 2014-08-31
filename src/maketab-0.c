/* ============================= C MeatAxe ==================================
   File:        $Id: maketab-0.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Calculate arithmetic tables
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include <string.h>
#include <stdlib.h>



#define MAXGRAD 12		/* Maximal degree of polynomials */
#define MAXSUBFIELDORD 16	/* Maximal order of subfields */
#define MAXSUBFIELDS 4		/* Maximal number of subfields */

typedef unsigned char BYTE;
typedef unsigned char POLY[MAXGRAD+1];


/* -----------------------------------------------------------------
   Global data
   ----------------------------------------------------------------- */

MTX_DEFINE_FILE_INFO

static MtxApplicationInfo_t AppInfo =
{
"maketab", "Generate table files (small fields version)",
"SYNTAX\n"
"    maketab " MTX_COMMON_OPTIONS_SYNTAX " <Field>\n"
"\n"
"ARGUMENTS\n"
"    <Field> ................. Field order (2..256)\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"\n"
"FILES\n"
"    p<Field>.zzz ............ O Table file\n"
};

static MtxApplication_t *App = NULL;

BYTE	tmult[256][256],
	tadd[256][256],
	tffirst[256][2],
	textract[8][256],
	taddinv[256],
	tmultinv[256],
	tnull[8][256],
	tinsert[8][256];
BYTE embed[MAXSUBFIELDS][MAXSUBFIELDORD]; /* Embeddings of subfields */
BYTE ffrestrict[MAXSUBFIELDS][256];	  /* Restriction to subfields */
long embedord[MAXSUBFIELDS];		  /* Subfield orders */

long info[4] = {0L,0L,0L,0L};
long ver = ZZZVERSION;
char filename[50];

long 	P;		/* Characteristic of the field */
long	G;		/* Generator for the field */
long	Q;		/* Order of the field */
long	CPM;		/* No. of field elements (FELs) per BYTE */
long	N;		/* Q = P^N */
long	maxmem;		/* (Highest value stored in BYTE) + 1 */
FILE	*fd;		/* File pointer */


POLY irred;		/*  Polynomial which defines the field */
BYTE indx[256];		/*  Index i of a field element g,g=X^i */
BYTE polynom[256];	/*  reverse to index  */
BYTE zech[256];		/*  Zech-logarithm for index i  */


/* Tables for non-prime fields, q<=256
   ----------------------------------- */

POLY irreducibles[] = {			/* Parker's polynomials: */
	{0,0,0,0,0,0,0,0,0,0,1,1,1},    /* F4   X2+X+1        */
	{0,0,0,0,0,0,0,0,0,1,0,1,1},    /* F8   X3+X+1        */
	{0,0,0,0,0,0,0,0,0,0,1,2,2},    /* F9   X2+2X+2       */
	{0,0,0,0,0,0,0,0,1,0,0,1,1},    /* F16  X4+X+1        */
	{0,0,0,0,0,0,0,0,0,0,1,4,2},    /* F25  X2+4X+2       */
	{0,0,0,0,0,0,0,0,0,1,0,2,1},    /* F27  X3+2X+1       */
	{0,0,0,0,0,0,0,1,0,0,1,0,1},    /* F32  X5+X2+1       */
	{0,0,0,0,0,0,0,0,0,0,1,6,3},    /* F49  X2+6X+3       */
	{0,0,0,0,0,0,1,0,1,1,0,1,1},    /* F64  X6+X4+X3+X+1  */
	{0,0,0,0,0,0,0,0,1,2,0,0,2},    /* F81  X4+2X3+2      */
	{0,0,0,0,0,0,0,0,0,0,1,7,2},    /* F121 X2+7X+2       */
	{0,0,0,0,0,0,0,0,0,1,0,3,3},    /* F125 X3+3X+3       */
	{0,0,0,0,0,1,0,0,0,0,0,1,1},    /* F128 X7+X+1        */
	{0,0,0,0,0,0,0,0,0,0,1,12,2},   /* F169 X2+12X+2      */
	{0,0,0,0,0,0,0,1,0,0,0,2,1},    /* F243 X5+2X+1       */
	{0,0,0,0,1,0,0,0,1,1,1,0,1}     /* F256 X8+X4+X3+X2+1 */
	};

/* The following tables contain the corresponding field orders
   and prime field orders for each of the polynomials above */

int irrednrs[] =	/* Field orders */
	{4,8,9,16,25,27,32,49,64,81,121,125,128,169,243,256,0};

BYTE irredprs[] =	/* Prime field orders  */
	{2,2,3,2,5,3,2,7,2,3,11,5,2,13,3,2,0};

/* The following is a list of possible generators for PRIME
   fields. For non-prime fields, X will be used as generator */

 BYTE gen[] = {1,2,3,5,6,7,19,0};


/* -----------------------------------------------------------------
   printpol() - Print a polynomial
   ----------------------------------------------------------------- */

void printpol(POLY a)

{
    int i,flag = 0;

    for (i = MAXGRAD; i >= 0; i--)
    {
	if (a[i]!=0)
	{	if (flag) printf("+");
		if (a[i] != 1) printf("%d",(int)a[i]);
		printf("x^%d",i);
		flag=1;
       	}
    }
    printf("\n");
}


/* -------------------------------------------------------------------------
   number() - Convert polynomial to number

   Given a polynomial f(x) over GF(p), this function interprets f(x) as a
   polynomial over Z, and returns the integer f(p).
   -------------------------------------------------------------------------- */

BYTE number(POLY a)

{
    BYTE k;
    int i;

    k = 0;
    for (i = MAXGRAD; i >= 0; i--)
	k = (BYTE) (k * P + a[i]);
    return k;
}


/* -----------------------------------------------------------------
   polmultx() - Multiply a polynomial by X.
   ----------------------------------------------------------------- */

void polmultx(POLY a)

{	int i;

	for (i = MAXGRAD; i > 0; --i)
		a[i] = a[i-1];
	a[0] = 0;
}


/* -----------------------------------------------------------------
   polymod() - Reduce the polynomial a modulo b. b ist assumed to
	be normalized.
   ----------------------------------------------------------------- */

static void polymod(POLY a, POLY b)

{
    int i, l, dl, f;

    /* l= index of leading coeff. in b (must be 1) */
    for (l = MAXGRAD; b[l]==0 && l>0; l--);
    for (dl = MAXGRAD; dl>=l; dl--)
    {	f = (int) a[dl];
	if (f == 0) continue;
	f = (int)P - f;
	for (i = 0; i <= l; ++i)
	    a[i+dl-l] = (BYTE) ((f*b[i] + a[i+dl-l]) % (int)P);
    }
}


/* -----------------------------------------------------------------
   testprim() - Test for primitivity.
   ----------------------------------------------------------------- */

void testprim()

{
    int i, a[256];

    memset(a,0,sizeof(a));
    for (i = 0; i < (int) Q; i++)
	a[indx[i]] += 1;
    for (i = 0; i < (int) Q; i++)
       	if(a[i] != 1)
	{
	    fprintf(stderr,"*** a[%d]=%d.",i,a[i]);
	    MTX_ERROR("Polynome is not primitive.");
	}
}


/* -----------------------------------------------------------------
   initarith() - Initialize index and zech logarithm tables.
   ----------------------------------------------------------------- */

void initarith()

{	int i,elem;
	POLY a;

	memset(indx,0,sizeof(indx));
	memset(a,0,sizeof(POLY));

	/* Initialize index table
	   ---------------------- */
	indx[0] = (BYTE) (Q - 1);
	polynom[(int)Q-1] = 0;		/* 0 gets index q-1 */
	a[0] = 1;			/* a=X^0  */
	for (i = 0; i < (int)Q-1; i++)	/* for each index */
	{	elem = number(a);
		indx[elem] = (BYTE) i;
		polynom[i] = (BYTE) elem;
		polmultx(a);
		polymod(a,irred);
        }
	testprim();

	/* Calculate zech logarithms
	   ------------------------- */
	for (i = 0; i <= (int)Q-1; i++)	/* for all elements (0 too) */
	{	elem = (int)((i%P)==P-1 ? i+1-P : i+1); /* add 1 */
		zech[indx[i]]=indx[elem]; /* Zech-table=result */
        }
}


/* -----------------------------------------------------------------
   add() - Add two field elements.
   ----------------------------------------------------------------- */

BYTE add(BYTE i, BYTE j)

{
    int ii,ij,z;

    if (P==Q) return ((BYTE) ( ((int)i+(int)j) % P));
    if (i==0) return(j);
    if (j==0) return(i);

    ii = indx[i];
    ij = indx[j];      /* x^a+x^b=x^(a+zech[(b-a)mod (q-1)])mod (q-1) */
    z = zech[(ij-ii+(int)Q-1) % ((int)Q-1)];
    if (z == (int)Q-1)
	return(0);             /* Zech-logarithm 0 */
    return (polynom[(ii+z) % ((int)Q-1)]);
}


/* -----------------------------------------------------------------
   mult() - Multiply two field elements.
   ----------------------------------------------------------------- */

BYTE mult(BYTE i, BYTE j)

{	if (P==Q)
		return ((BYTE)(((long int)i * j) % P));
	if (i==0 || j==0)
		return(0);

	return (polynom[(indx[i] + indx[j]) % ((int)Q-1)]);
}


/* ------------------------------------------------------------------
   testgen() - Test if ord(a) = prime-1
   ------------------------------------------------------------------ */

int testgen(BYTE a, BYTE prime)

{	BYTE i, x;

	if (a % prime == 0) return 0;
	x = a;
	for (i = 1; x != 1; ++i)
	{	x = (BYTE) (((long int)x * a) % prime);
	}
	return (i == prime - (BYTE) 1);
}


/* ------------------------------------------------------------------
   unpack() - Unpack a BYTE into an array of field elements
   pack() - Pack field elements into one BYTE

   We use q-adic packing, i.e.

	pack(a[0]...a[n]) := a[n]*q^0 + ... + a[0]*q^n

   with n = CPM - 1.
   ------------------------------------------------------------------ */

void unpack(register BYTE x, BYTE a[8])

{	int i;

	for (i = (int)(CPM-1); i >= 0; i--)
	{	a[i] = (BYTE) ((int) x % (int) Q);
		x = (BYTE)((int) x / (int) Q);
	}
}


BYTE pack(BYTE a[8])

{	int i;
	BYTE x;

	x = 0;
	for (i = 0; i < (int) CPM; i++)
		x = (BYTE)(x * Q + a[i]);
	return (x);
}


/* -----------------------------------------------------------------
   writeheader() - Set info[], open table file, select polynomial,
	and initialize tables.
   ----------------------------------------------------------------- */

void writeheader()

{
    int i, j;

    sprintf(filename,"p%3.3ld.zzz",Q);
    fd = SysFopen(filename,FM_CREATE);
    if (fd == NULL)
    {
	perror(filename);
	MTX_ERROR("Cannot open table file");
    }
    for (CPM=1,maxmem=Q; (long)maxmem * Q <= 256L; ++CPM, maxmem *= Q);
    for (i = 0; irrednrs[i] != (int) Q && irrednrs[i] != 0; ++i);
    if (irrednrs[i] != 0)
    {
	P = irredprs[i];
        for (j = 0; j <= MAXGRAD; j++)
            irred[j] = irreducibles[i][MAXGRAD-j];
	G = P;		/* Generator is X */
	initarith();	/* Init index- and Zech-tables */
    }
    else
    {	
	P = Q;
	/* Find a generator
	   ---------------- */
	for (j = 0; (G = (long) gen[j]) != 0 &&
		 !testgen((BYTE)gen[j],(BYTE)P); j++);
    }
    info[0] = (long int) P;
    info[1] = (long int) G;
    info[2] = (long int) Q;
    info[3] = (long int) CPM;

    MESSAGE(1,("ZZZ version : %ld\n",ver));
    MESSAGE(1,("Field order : %ld=%ld^%ld\n",info[2],info[0],N));
    if (P != Q && MtxMessageLevel >= 1)
    {
	printf("Polynome    : ");
	printpol(irred);
    }
    MESSAGE(1,("Generator   : %ld\n",info[1]));
    MESSAGE(1,("Packing     : %ld/byte\n",info[3]));
}


/* -----------------------------------------------------------------
   checkq() - Set Q and N. Verify that Q is a prime power.
   ----------------------------------------------------------------- */

void checkq(long l)

{
    long q, d;

    if (l < 2 || l > 256)
    {
	fprintf(stderr,"Field order out of range (2-256)\n");
	exit(EXIT_ERR);
    }

    Q = l;
    q = Q;
    for (d = 2; q % d != 0; ++d); /* Find smallest prime divisor */
    for (N = 0; (q % d) == 0; ++N)
       	q /= d;
    if (q != 1)
    {
	fprintf(stderr,"Illegal Field order\n");
	exit(EXIT_ERR);
    }
}


/* -----------------------------------------------------------------
   inittables() - Initialize arithmetic tables with 0xFF
   ----------------------------------------------------------------- */

void inittables()

{
	memset(tmult,0xFF,sizeof(tmult));
	memset(tadd,0xFF,sizeof(tadd));
	memset(tffirst,0xFF,sizeof(tffirst));
	memset(textract,0xFF,sizeof(textract));
	memset(taddinv,0xFF,sizeof(taddinv));
	memset(tmultinv,0xFF,sizeof(tmultinv));
	memset(tnull,0xFF,sizeof(tnull));
	memset(tinsert,0xFF,sizeof(tinsert));
}

/* -----------------------------------------------------------------
   mkembed() - Calculate embeddings of all subfields.
   ----------------------------------------------------------------- */

void mkembed()

{
    int n;	/* Degree of subfield over Z_p */
    long q; /* subfield order */
    int i, k;
    POLY a, subirred;
    int count = 0;
    BYTE emb, f;

    memset(embed,255,sizeof(embed));
    memset(ffrestrict,255,sizeof(ffrestrict));

    MESSAGE(1,("Calculating embeddings of subfields\n"));

    /* Clear the embedord array. embedord[i]=0 means
       that the entry (and all subequent) is not used.
       ----------------------------------------------- */
    for (i = 0; i < MAXSUBFIELDS; ++i) embedord[i] = 0;


    for (n = 1; n < N; ++n)
    {
	if (N % n != 0) continue;	/* n must divide N */

        /* The prime field is simple:
           -------------------------- */
	if (n == 1)
	{
	    MESSAGE(1,("GF(%ld)\n",P));
    	    embedord[count] = P;
    	    for (i = 0; i < (int) P; ++i)
    	    {
		embed[count][i] = (BYTE) i;
		ffrestrict[count][i] = (BYTE) i;
    	    }
	    ++count;
	    continue;
	}

	/* Calculate the subfield order
	   ---------------------------- */
	for (q = 1, i = n; i > 0; --i, q *= P);
	embedord[count] = q;
	embed[count][0] = 0;
	ffrestrict[count][0] = 0;
	if ((Q-1) % (q-1) != 0)
	{
	    fprintf(stderr,"*** q=%ld, Q=%ld.",q,Q);
	    MTX_ERROR("Internal error.");
	}

	/* Calculate a generator for the subfield
	   -------------------------------------- */
	emb = 1;
	for (i = (Q-1)/(q-1); i > 0; --i) 
	    emb = mult(emb,(BYTE) G);

	/* Look up the polynomial
	   ---------------------- */
	for (k = 0; irrednrs[k] != 0 && irrednrs[k] != q; ++k);
	if (irrednrs[k] == 0) MTX_ERROR("Internal error.");
       	for (i = 0; i <= MAXGRAD; i++)
            subirred[i] = irreducibles[k][MAXGRAD-i];

	MESSAGE(1,("GF(%ld): gen=%d pol=",q,emb));
	if (MSG1) printpol(subirred);
	fflush(stdout);

	memset(a,0,sizeof(POLY));
	a[0] = 1;		/* a=X^0  */
	f = FF_ONE;
	for (i = 0; i < (int)q-1; ++i)
	{
	    embed[count][number(a)] = f;
	    MESSAGE(3,("embed[%d][%d]=%d\n",count,number(a),(int)f));
	    ffrestrict[count][f] = number(a);
	    polmultx(a);
	    polymod(a,subirred);
	    f = mult(f,emb);
        }
	++count;
    }
    MESSAGE(1,("\n"));

    if (MtxMessageLevel >= 2)
    {
        for (i = 0; i < 4; ++i)
        {
	    printf("  GF(%2ld): ",embedord[i]);
            for (k=0; k < 16; ++k)
	        printf("%4d",embed[i][k]);
 	    printf("\n");
	    fflush(stdout);
	}
    }
}


static int Init(int argc, const char **argv)

{
    long l;
    char rev[10], date[20];

    if ((App = AppAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    if (AppGetArguments(App,1,1) < 0)
	return -1;
    if (sscanf(App->ArgV[0],"%ld",&l) != 1)
	MTX_ERROR2("%s: %E\n",AppInfo.Name,MTX_ERR_BADUSAGE);
    checkq(l);
    sscanf(MtxVersion,"%s %*d %s",rev,date);
    MESSAGE(0,("MAKETAB Revision %s (%s)\n",rev,date));
    return 0;
}

/* -----------------------------------------------------------------
   main()
   ----------------------------------------------------------------- */

int main(int argc, char *argv[])

{
    int i, j, k;
    short flag;
    BYTE a[8],b[8],c[8],d[8],z;


    /* Initialize
       ---------- */
    if (Init(argc,(const char **)argv) != 0)
	return 1;
    writeheader();			/* Open file and write header */
    inittables();

    /* Make insert table
       ----------------- */
    memset(a,0,sizeof(a));
    MESSAGE(1,("Calculating insert table\n"));
    for (i = 0; i < (int) Q; i++)
    {
	for (j = 0; j < (int) CPM; j++)
    	{
	    a[j] = (BYTE) i;
	    tinsert[j][i] = pack(a);	/* Insert-table */
	    MESSAGE(3,("insert[%d][%d]=%u (0x%x)\n",j,i,
	     tinsert[j][i],tinsert[j][i]));
	    a[j] = 0;
    	}
    }

    /* Pack/unpack and arithmetic tables
       --------------------------------- */
    MESSAGE(1,("Calculating pack/unpack and arithmetic tables\n"));
    for (i=0 ; i < (int) maxmem; i++)
    {
	if (i % 10 == 0 && MtxMessageLevel >= 2)
	{   if (i == 140) printf("\n");
	    printf("%3d ",i);
	}
	unpack((BYTE)i,a);   /* unpack 1.element in a[] */
	flag = 0;
	for (j = 0; j < (int) CPM; j++)
	{
	    textract[j][i] = a[j];
	    z = a[j];
	    a[j] = 0;
	    tnull[j][i] = pack(a);     /* Null-table */
	    a[j] = z;
	    if (!flag && z)
	    {
		flag = 1;
		tffirst[i][0] = z;  /* Find first table: mark */
		tffirst[i][1] = (BYTE)j;  /* Find first table: pos. */
	    }
	}
	if (Q != 2)
	{
	    for (j=0; j < (int) maxmem; j++)
	    {
		unpack((BYTE)j,b);	/* 2.element in b[] */
		if (i <= j)
		{
		    for (k=0; k < (int) CPM; k++)
			c[k] = add(a[k],b[k]);
		    tadd[i][j] = pack(c);
		}
		else
		    tadd[i][j]=tadd[j][i];

		if (i < (int) Q)
		{
		    for (k=0; k < (int) CPM; k++)
			d[k] = mult(a[(int)CPM-1],b[k]);
		    tmult[i][j] = pack(d);
		}
		else
		    tmult[i][j] = tmult[i-(int)Q][j];
	    }
	}
	else	/* GF(2) */
	{
	    for (j=0; j < (int) maxmem; j++)
	    {
		tadd[i][j] = (BYTE)(i ^ j);
		tmult[i][j] = (BYTE)((i & 1) != 0 ?  j : 0);
	    }
	}
    }
    MESSAGE(2,("\n"));

    /* Inversion table
       --------------- */
    MESSAGE(1,("Calculating inversion table\n"));
    fflush(stdout);
    for (i = 0; i < (int)Q; i++)
    {
        for (j = 0; j < (int)Q; j++)
	{
	    if (add((BYTE)i,(BYTE)j) == 0) taddinv[i] = (BYTE)j;
	    if (mult((BYTE)i,(BYTE)j) == 1) tmultinv[i] = (BYTE)j;
	}
    }

    mkembed();

    MESSAGE(0,("Writing tables to %s\n",filename));
    if (
      SysWriteLong(fd,info,4) != 4 ||
      SysWriteLong(fd,&ver,1) != 1 ||
      fwrite(tmult,4,0x4000,fd) != 0x4000 ||
      fwrite(tadd,4,0x4000,fd) != 0x4000 ||
      fwrite(tffirst,1,sizeof(tffirst),fd) != sizeof(tffirst) ||
      fwrite(textract,1,sizeof(textract),fd)!= sizeof(textract) ||
      fwrite(taddinv,1,sizeof(taddinv),fd) != sizeof(taddinv) ||
      fwrite(tmultinv,1,sizeof(tmultinv),fd)!= sizeof(tmultinv) ||
      fwrite(tnull,1,sizeof(tnull),fd) != sizeof(tnull) ||
      fwrite(tinsert,1,sizeof(tinsert),fd) != sizeof(tinsert) ||
      SysWriteLong(fd,embedord,MAXSUBFIELDS) != MAXSUBFIELDS ||
      fwrite(embed,MAXSUBFIELDORD,MAXSUBFIELDS,fd) != MAXSUBFIELDS ||
      fwrite(ffrestrict,256,MAXSUBFIELDS,fd) != MAXSUBFIELDS
      )
    {
	perror(filename);
	MTX_ERROR("Error writing table file");
    }
    fclose(fd);
    AppFree(App);
    return(0);
}



/**
@page prog_maketab maketab - Calculate Arithmetic Tables

<<<<<<< HEAD
@section syntax Command Line
=======
@section maketab-0_syntax Command Line
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
<pre>
maketab @em Field
</pre>

@par @em Field
  The field order.

<<<<<<< HEAD
@section out Output Files
@par pXXX.zzz
  Arithmetic tables for GF(XXX).

@section desc Description
=======
@section maketab-0_out Output Files
@par pXXX.zzz
  Arithmetic tables for GF(XXX).

@section maketab-0_desc Description
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c

@attention
   This program is no longer needed. Arithmetic tables are calculated automatically
   when needed.


This program generates the lookup tables used by the finite field
arithmetic. Tables are stored in a file (one file for each field) 
which is loaded by the MeatAxe applications at run-time. The argument,
@em Field, must be a positive integer (a prime power, actually). The 
upper limit for the field order depends on which version of the arithmetic 
you are using. The arithmetic version is chosen at compile-time.

The standard version of MAKETAB creates table files for fields up to 
GF(256). The output file is named "pXXX.zzz", where XXX is the field order.
The size of this file is approximately 130 KBytes and does not depend on
the field order. For the `big' version, the largest possible field is 
GF(65536), and the output file is named "pXXXXX.zzz", where XXXXX is 
the field order. This file is usually (i.e., for comparable field orders) 
much smaller than the standard version file, but its size increases 
linearly with the field order.

The table file is needed whenever a program performs finite field
operations, i.e., when working with matrices. The programs look for
the table file first in the current directory and then, if it is not
there, in the library directory. Thus,
table files for the most commonly used fields can be kept in the
library, and additional table files can be created in the current
directory when they are needed.
Usually, the name of the library directory is defined at compile-time,
but this may be changed at any time by defining the environment
variable MTXLIB. 

*/

