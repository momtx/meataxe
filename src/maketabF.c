////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Calculate arithmetic tables
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////


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

BYTE
    mtx_tmult[256][256],
    mtx_tadd[256][256],
    mtx_taddinv[256],
    mtx_tmultinv[256],
	mtx_tffirst[256][2],
	mtx_textract[8][256],
	mtx_tnull[8][256],
	mtx_tinsert[8][256];
BYTE mtx_embed[MAXSUBFIELDS][MAXSUBFIELDORD]; /* Embeddings of subfields */
BYTE mtx_restrict[MAXSUBFIELDS][256];	  /* Restriction to subfields */
long mtx_embedord[MAXSUBFIELDS];		  /* Subfield orders */

static long info[4] = {0L,0L,0L,0L};
static long ver = ZZZVERSION;
static char filename[50];

static long 	P;		/* Characteristic of the field */
static long	G;		/* Generator for the field */
static long	Q;		/* Order of the field */
static long	CPM;		/* No. of field elements (FELs) per BYTE */
static long	N;		/* Q = P^N */
static long	maxmem;		/* (Highest value stored in BYTE) + 1 */
static FILE	*fd;		/* File pointer */


static POLY irred;		/*  Polynomial which defines the field */
static BYTE indx[256];		/*  Index i of a field element g,g=X^i */
static BYTE polynom[256];	/*  reverse to index  */
static BYTE zech[256];		/*  Zech-logarithm for index i  */


/* Tables for non-prime fields, q<=256
   ----------------------------------- */

static POLY irreducibles[] = {			/* Parker's polynomials: */
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

static int irrednrs[] =	/* Field orders */
	{4,8,9,16,25,27,32,49,64,81,121,125,128,169,243,256,0};

static BYTE irredprs[] =	/* Prime field orders  */
	{2,2,3,2,5,3,2,7,2,3,11,5,2,13,3,2,0};

/* The following is a list of possible generators for PRIME
   fields. For non-prime fields, X will be used as generator */

static BYTE gen[] = {1,2,3,5,6,7,19,0};


/* -----------------------------------------------------------------
   printpol() - Print a polynomial
   ----------------------------------------------------------------- */

static void printpol(POLY a)
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

static BYTE number(POLY a)
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

static void polmultx(POLY a)
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

static void testprim()
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

static void initarith()
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

static BYTE add(BYTE i, BYTE j)
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

static BYTE mult(BYTE i, BYTE j)

{	if (P==Q)
		return ((BYTE)(((long int)i * j) % P));
	if (i==0 || j==0)
		return(0);

	return (polynom[(indx[i] + indx[j]) % ((int)Q-1)]);
}


/* ------------------------------------------------------------------
   testgen() - Test if ord(a) = prime-1
   ------------------------------------------------------------------ */

static int testgen(BYTE a, BYTE prime)
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

static void unpack(register BYTE x, BYTE a[8])
{	int i;

	for (i = (int)(CPM-1); i >= 0; i--)
	{	a[i] = (BYTE) ((int) x % (int) Q);
		x = (BYTE)((int) x / (int) Q);
	}
}


static BYTE pack(BYTE a[8])
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

static void writeheader()
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

static void checkq(long l)
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

static void inittables()
{
	memset(mtx_tmult,0xFF,sizeof(mtx_tmult));
	memset(mtx_tadd,0xFF,sizeof(mtx_tadd));
	memset(mtx_tffirst,0xFF,sizeof(mtx_tffirst));
	memset(mtx_textract,0xFF,sizeof(mtx_textract));
	memset(mtx_taddinv,0xFF,sizeof(mtx_taddinv));
	memset(mtx_tmultinv,0xFF,sizeof(mtx_tmultinv));
	memset(mtx_tnull,0xFF,sizeof(mtx_tnull));
	memset(mtx_tinsert,0xFF,sizeof(mtx_tinsert));
}

/* -----------------------------------------------------------------
   mkembed() - Calculate embeddings of all subfields.
   ----------------------------------------------------------------- */

static void mkembed()
{
    int n;	/* Degree of subfield over Z_p */
    long q; /* subfield order */
    int i, k;
    POLY a, subirred;
    int count = 0;
    BYTE emb, f;

    memset(mtx_embed,255,sizeof(mtx_embed));
    memset(mtx_restrict,255,sizeof(mtx_restrict));

    MESSAGE(1,("Calculating embeddings of subfields\n"));

    /* Clear the mtx_embedord array. mtx_embedord[i]=0 means
       that the entry (and all subequent) is not used.
       ----------------------------------------------- */
    for (i = 0; i < MAXSUBFIELDS; ++i) mtx_embedord[i] = 0;


    for (n = 1; n < N; ++n)
    {
	if (N % n != 0) continue;	/* n must divide N */

        /* The prime field is simple:
           -------------------------- */
	if (n == 1)
	{
	    MESSAGE(1,("GF(%ld)\n",P));
    	    mtx_embedord[count] = P;
    	    for (i = 0; i < (int) P; ++i)
    	    {
		mtx_embed[count][i] = (BYTE) i;
		mtx_restrict[count][i] = (BYTE) i;
    	    }
	    ++count;
	    continue;
	}

	/* Calculate the subfield order
	   ---------------------------- */
	for (q = 1, i = n; i > 0; --i, q *= P);
	mtx_embedord[count] = q;
	mtx_embed[count][0] = 0;
	mtx_restrict[count][0] = 0;
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
	    mtx_embed[count][number(a)] = f;
	    MESSAGE(3,("mtx_embed[%d][%d]=%d\n",count,number(a),(int)f));
	    mtx_restrict[count][f] = number(a);
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
	    printf("  GF(%2ld): ",mtx_embedord[i]);
            for (k=0; k < 16; ++k)
	        printf("%4d",mtx_embed[i][k]);
 	    printf("\n");
	    fflush(stdout);
	}
    }
}


static int Init(int field)
{
    checkq(field);
    return 0;
}

/* -----------------------------------------------------------------
   FfMakeTables() - Build meataxe arithmetic tables
   ----------------------------------------------------------------- */

int FfMakeTables(int field)
{
    int i, j, k;
    short flag;
    BYTE a[8],b[8],c[8],d[8],z;


    /* Initialize
       ---------- */
    if (Init(field) != 0)
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
	    mtx_tinsert[j][i] = pack(a);	/* Insert-table */
	    MESSAGE(3,("insert[%d][%d]=%u (0x%x)\n",j,i,
	     mtx_tinsert[j][i],mtx_tinsert[j][i]));
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
	    mtx_textract[j][i] = a[j];
	    z = a[j];
	    a[j] = 0;
	    mtx_tnull[j][i] = pack(a);     /* Null-table */
	    a[j] = z;
	    if (!flag && z)
	    {
		flag = 1;
		mtx_tffirst[i][0] = z;  /* Find first table: mark */
		mtx_tffirst[i][1] = (BYTE)j;  /* Find first table: pos. */
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
		    mtx_tadd[i][j] = pack(c);
		}
		else
		    mtx_tadd[i][j]=mtx_tadd[j][i];

		if (i < (int) Q)
		{
		    for (k=0; k < (int) CPM; k++)
			d[k] = mult(a[(int)CPM-1],b[k]);
		    mtx_tmult[i][j] = pack(d);
		}
		else
		    mtx_tmult[i][j] = mtx_tmult[i-(int)Q][j];
	    }
	}
	else	/* GF(2) */
	{
	    for (j=0; j < (int) maxmem; j++)
	    {
		mtx_tadd[i][j] = (BYTE)(i ^ j);
		mtx_tmult[i][j] = (BYTE)((i & 1) != 0 ?  j : 0);
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
	    if (add((BYTE)i,(BYTE)j) == 0) mtx_taddinv[i] = (BYTE)j;
	    if (mult((BYTE)i,(BYTE)j) == 1) mtx_tmultinv[i] = (BYTE)j;
	}
    }

    mkembed();

    MESSAGE(1,("Writing tables to %s\n",filename));
    if (
      SysWriteLong(fd,info,4) != 4 ||
      SysWriteLong(fd,&ver,1) != 1 ||
      fwrite(mtx_tmult,4,0x4000,fd) != 0x4000 ||
      fwrite(mtx_tadd,4,0x4000,fd) != 0x4000 ||
      fwrite(mtx_tffirst,1,sizeof(mtx_tffirst),fd) != sizeof(mtx_tffirst) ||
      fwrite(mtx_textract,1,sizeof(mtx_textract),fd)!= sizeof(mtx_textract) ||
      fwrite(mtx_taddinv,1,sizeof(mtx_taddinv),fd) != sizeof(mtx_taddinv) ||
      fwrite(mtx_tmultinv,1,sizeof(mtx_tmultinv),fd)!= sizeof(mtx_tmultinv) ||
      fwrite(mtx_tnull,1,sizeof(mtx_tnull),fd) != sizeof(mtx_tnull) ||
      fwrite(mtx_tinsert,1,sizeof(mtx_tinsert),fd) != sizeof(mtx_tinsert) ||
      SysWriteLong(fd,mtx_embedord,MAXSUBFIELDS) != MAXSUBFIELDS ||
      fwrite(mtx_embed,MAXSUBFIELDORD,MAXSUBFIELDS,fd) != MAXSUBFIELDS ||
      fwrite(mtx_restrict,256,MAXSUBFIELDS,fd) != MAXSUBFIELDS
      )
    {
	perror(filename);
	MTX_ERROR("Error writing table file");
    }
    fclose(fd);
    return(0);
}

