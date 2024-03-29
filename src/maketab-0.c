////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Calculate arithmetic tables
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"
#include <string.h>
#include <stdlib.h>



#define MAXGRAD 12		/* Maximal degree of polynomials */

typedef unsigned char POLY[MAXGRAD+1];


/* -----------------------------------------------------------------
   Global data
   ----------------------------------------------------------------- */


uint8_t
    mtx_tmult[256][256],
    mtx_tadd[256][256],
    mtx_taddinv[256],
    mtx_tmultinv[256],
	mtx_tffirst[256][2],
	mtx_textract[8][256],
	mtx_tnull[8][256],
	mtx_tinsert[8][256];
uint8_t mtx_embed[MTX_MAXSUBFIELDS][MTX_MAXSUBFIELDORD]; /* Embeddings of subfields */
uint8_t mtx_restrict[MTX_MAXSUBFIELDS][256];	  /* Restriction to subfields */
static uint32_t subfieldOrder[MTX_MAXSUBFIELDS];		  /* Subfield orders */

static uint32_t info[4] = {0L,0L,0L,0L};
static uint32_t ver = MTX_ZZZVERSION;
static char filename[50];

static long Q;		// field order
static long P;		// characteristic
static long N;		// degree over prime field, Q = P^N
static long G;		// generator for the multiplicative group
static long CPM;	// no. of field elements (FELs) per uint8_t
static long maxmem;	// (highest value stored in uint8_t) + 1
static FILE *fd;	// table file pointer


static POLY irred;		/*  Polynomial which defines the field */
static uint8_t indx[256];		/*  Index i of a field element g,g=X^i */
static uint8_t polynom[256];	/*  reverse to index  */
static uint8_t zech[256];		/*  Zech-logarithm for index i  */


/* Tables for non-prime fields, q<=256
   ----------------------------------- */

static POLY irreducibles[] = {		// Parker's polynomials:
	{0,0,0,0,0,0,0,0,0,0,1,1,1},    // F4   X2+X+1
	{0,0,0,0,0,0,0,0,0,1,0,1,1},    // F8   X3+X+1
	{0,0,0,0,0,0,0,0,0,0,1,2,2},    // F9   X2+2X+2
	{0,0,0,0,0,0,0,0,1,0,0,1,1},    // F16  X4+X+1
	{0,0,0,0,0,0,0,0,0,0,1,4,2},    // F25  X2+4X+2
	{0,0,0,0,0,0,0,0,0,1,0,2,1},    // F27  X3+2X+1
	{0,0,0,0,0,0,0,1,0,0,1,0,1},    // F32  X5+X2+1
	{0,0,0,0,0,0,0,0,0,0,1,6,3},    // F49  X2+6X+3
	{0,0,0,0,0,0,1,0,1,1,0,1,1},    // F64  X6+X4+X3+X+1
	{0,0,0,0,0,0,0,0,1,2,0,0,2},    // F81  X4+2X3+2
	{0,0,0,0,0,0,0,0,0,0,1,7,2},    // F121 X2+7X+2 
	{0,0,0,0,0,0,0,0,0,1,0,3,3},    // F125 X3+3X+3  
	{0,0,0,0,0,1,0,0,0,0,0,1,1},    // F128 X7+X+1    
	{0,0,0,0,0,0,0,0,0,0,1,12,2},   // F169 X2+12X+2   
	{0,0,0,0,0,0,0,1,0,0,0,2,1},    // F243 X5+2X+1     
	{0,0,0,0,1,0,0,0,1,1,1,0,1}     // F256 X8+X4+X3+X2+1
	};

/* The following tables contain the corresponding field orders
   and prime field orders for each of the polynomials above */

static int irrednrs[] =	/* Field orders */
	{4,8,9,16,25,27,32,49,64,81,121,125,128,169,243,256,0};

static uint8_t irredprs[] =	/* Prime field orders  */
	{2,2,3,2,5,3,2,7,2,3,11,5,2,13,3,2,0};

/* The following is a list of possible generators for PRIME
   fields. For non-prime fields, X will be used as generator */

static uint8_t gen[] = {1,2,3,5,6,7,19,0};



////////////////////////////////////////////////////////////////////////////////////////////////////

static void formatPoly(StrBuffer_t* buf, POLY a)
{
   int i, flag = 0;

   for (i = MAXGRAD; i >= 0; i--) {
      if (a[i] != 0) {
         if (flag) { sbAppend(buf, "+"); }
         if (a[i] != 1) { sbPrintf(buf, "%d", (int)a[i]); }
         sbPrintf(buf, "x^%d", i);
         flag = 1;
      }
   }
}

/* -------------------------------------------------------------------------
   number() - Convert polynomial to number

   Given a polynomial f(x) over GF(p), this function interprets f(x) as a
   polynomial over Z, and returns the integer f(p).
   -------------------------------------------------------------------------- */

static uint8_t number(POLY a)
{
    uint8_t k;
    int i;

    k = 0;
    for (i = MAXGRAD; i >= 0; i--)
	k = (uint8_t) (k * P + a[i]);
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
	    a[i+dl-l] = (uint8_t) ((f*b[i] + a[i+dl-l]) % (int)P);
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
	    mtxAbort(MTX_HERE,"Polynome is not primitive.");
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
	indx[0] = (uint8_t) (Q - 1);
	polynom[(int)Q-1] = 0;		/* 0 gets index q-1 */
	a[0] = 1;			/* a=X^0  */
	for (i = 0; i < (int)Q-1; i++)	/* for each index */
	{	elem = number(a);
		indx[elem] = (uint8_t) i;
		polynom[i] = (uint8_t) elem;
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

static uint8_t add(uint8_t i, uint8_t j)
{
    int ii,ij,z;

    if (P==Q) return ((uint8_t) ( ((int)i+(int)j) % P));
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

static uint8_t mult(uint8_t i, uint8_t j)

{	if (P==Q)
		return ((uint8_t)(((long int)i * j) % P));
	if (i==0 || j==0)
		return(0);

	return (polynom[(indx[i] + indx[j]) % ((int)Q-1)]);
}


/* ------------------------------------------------------------------
   testgen() - Test if ord(a) = prime-1
   ------------------------------------------------------------------ */

static int testgen(uint8_t a, uint8_t prime)
{	uint8_t i, x;

	if (a % prime == 0) return 0;
	x = a;
	for (i = 1; x != 1; ++i)
	{	x = (uint8_t) (((long int)x * a) % prime);
	}
	return (i == prime - (uint8_t) 1);
}


/* ------------------------------------------------------------------
   unpack() - Unpack a uint8_t into an array of field elements
   pack() - Pack field elements into one uint8_t

   We use q-adic packing, i.e.

	pack(a[0]...a[n]) := a[n]*q^0 + ... + a[0]*q^n

   with n = CPM - 1.
   ------------------------------------------------------------------ */

static void unpack(register uint8_t x, uint8_t a[8])
{	int i;

	for (i = (int)(CPM-1); i >= 0; i--)
	{	a[i] = (uint8_t) ((int) x % (int) Q);
		x = (uint8_t)((int) x / (int) Q);
	}
}


static uint8_t pack(uint8_t a[8])
{	int i;
	uint8_t x;

	x = 0;
	for (i = 0; i < (int) CPM; i++)
		x = (uint8_t)(x * Q + a[i]);
	return (x);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

///  Set info[], open table file, select polynomial, and initialize tables.

static void writeheader()
{
    int i, j;

    sprintf(filename,"p%3.3ld.zzz",Q);
    fd = sysFopen(filename,"wb::lib");
    if (fd == NULL)
    {
	perror(filename);
	mtxAbort(MTX_HERE,"Cannot open table file");
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
	for (j = 0; (G = (long) gen[j]) != 0 && !testgen((uint8_t)gen[j],(uint8_t)P); j++);
    }
    info[0] = (long int) P;
    info[1] = (long int) G;
    info[2] = (long int) Q;
    info[3] = (long int) CPM;

    MTX_LOGD("ZZZ version : %lu",(unsigned long)ver);
    MTX_LOGD("Field order : %lu=%lu^%lu",
             (unsigned long)info[2],(unsigned long)info[0],(unsigned long)N);
    if (P != Q) {
       MTX_XLOGD(msg) {
          sbAppend(msg, "Polynome    : ");
          formatPoly(msg, irred);
       }
    }
    MTX_LOGD("Generator   : %lu",(unsigned long)info[1]);
    MTX_LOGD("Packing     : %lu/byte",(unsigned long)info[3]);
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

////////////////////////////////////////////////////////////////////////////////////////////////////
   
// Calculate embeddings of all subfields.

static void mkembed()
{
   int n;   // degree of subfield over Z_p
   long q;  // subfield order
   int i, k;
   POLY a, subirred;
   int count = 0;  // number of proper subfields
   uint8_t emb, f;

   memset(mtx_embed, 255, sizeof(mtx_embed));
   memset(mtx_restrict, 255, sizeof(mtx_restrict));
   memset(subfieldOrder, 0, sizeof(subfieldOrder));       // mark as unused

   MTX_LOGD("Calculating embeddings of subfields");

   for (n = 1; n < N; ++n) {
      // All subfields of GF(p^N) have order are p^n with n|N.
      if (N % n != 0) { continue; }

      // The prime field embedding is simple
      if (n == 1) {
         MTX_LOGD("GF(%ld)", P);
         subfieldOrder[count] = P;
         for (i = 0; i < (int) P; ++i) {
            mtx_embed[count][i] = (uint8_t) i;
            mtx_restrict[count][i] = (uint8_t) i;
         }
         ++count;
         continue;
      }

      // Calculate the subfield order
      for (q = 1, i = n; i > 0; --i, q *= P) {}
      subfieldOrder[count] = q;
      mtx_embed[count][0] = 0;
      mtx_restrict[count][0] = 0;
      MTX_ASSERT((Q - 1) % (q - 1) == 0);

      // Calculate a generator for the subfield
      emb = 1;
      for (i = (Q - 1) / (q - 1); i > 0; --i) {
         emb = mult(emb, (uint8_t) G);
      }

      // Look up the polynomial
      for (k = 0; irrednrs[k] != 0 && irrednrs[k] != q; ++k) {}
      MTX_ASSERT(irrednrs[k] != 0);
      for (i = 0; i <= MAXGRAD; i++) {
         subirred[i] = irreducibles[k][MAXGRAD - i];
      }

      MTX_XLOGD(msg) {
         sbPrintf(msg, "GF(%ld): gen=%d pol=", q, emb);
         formatPoly(msg, subirred);
      }

      memset(a, 0, sizeof(POLY));
      a[0] = 1;                 // a=X^0
      f = FF_ONE;
      for (i = 0; i < (int)q - 1; ++i) {
         mtx_embed[count][number(a)] = f;
         mtx_restrict[count][f] = number(a);
         polmultx(a);
         polymod(a, subirred);
         f = mult(f, emb);
      }
      ++count;
   }

   for (i = 0; i < 4 && subfieldOrder[i] != 0; ++i) {
      MTX_XLOG2(msg) {
         sbPrintf(msg, "GF(%2d) embedding: ", (int)subfieldOrder[i]);
         for (k = 0; k < 16; ++k) {
            sbPrintf(msg, "%4d", mtx_embed[i][k]);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void Init(int field)
{
   if (field < 2 || field > 256) {
      mtxAbort(MTX_HERE, "Field order %d out of range (2-256)", field);
   }

   // Check that Q is a prime power
   Q = field;
   long q = Q;
   long p;
   for (p = 2; q % p != 0; ++p);
   for (N = 0; (q % p) == 0; ++N) {
      q /= p;
   }
   if (q != 1) {
      mtxAbort(MTX_HERE, "Illegal field order %d", field);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Create the arithmetic table file.

int ffMakeTables(int field)
{
    int i, j, k;
    uint16_t flag;
    uint8_t a[8],b[8],c[8],d[8],z;


    // Initialize
    Init(field);
    writeheader();
    inittables();

    // Make insert table
    memset(a,0,sizeof(a));
    MTX_LOGD("Calculating insert table");
    for (i = 0; i < (int) Q; i++)
    {
	for (j = 0; j < (int) CPM; j++)
    	{
	    a[j] = (uint8_t) i;
	    mtx_tinsert[j][i] = pack(a);	/* Insert-table */
	    MTX_LOG2("insert[%d][%d]=%u (0x%x)",j,i,
	     mtx_tinsert[j][i],mtx_tinsert[j][i]);
	    a[j] = 0;
    	}
    }

    // Pack/unpack and arithmetic tables
    MTX_LOGD("Calculating pack/unpack and arithmetic tables");
    for (i=0 ; i < (int) maxmem; i++)
    {
	unpack((uint8_t)i,a);   /* unpack 1.element in a[] */
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
		mtx_tffirst[i][1] = (uint8_t)j;  /* Find first table: pos. */
	    }
	}
	if (Q != 2)
	{
	    for (j=0; j < (int) maxmem; j++)
	    {
		unpack((uint8_t)j,b);	/* 2.element in b[] */
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
		mtx_tadd[i][j] = (uint8_t)(i ^ j);
		mtx_tmult[i][j] = (uint8_t)((i & 1) != 0 ?  j : 0);
	    }
	}
    }

    // Inversion table
    MTX_LOGD("Calculating inversion table");
    for (i = 0; i < (int)Q; i++)
    {
        for (j = 0; j < (int)Q; j++)
	{
	    if (add((uint8_t)i,(uint8_t)j) == 0) mtx_taddinv[i] = (uint8_t)j;
	    if (mult((uint8_t)i,(uint8_t)j) == 1) mtx_tmultinv[i] = (uint8_t)j;
	}
    }

    mkembed();

    MTX_LOGD("Writing tables to %s",filename);
    sysWrite32(fd,info,4);
    sysWrite32(fd,&ver,1);
    sysWrite8(fd,mtx_tmult,sizeof(mtx_tmult));
    sysWrite8(fd,mtx_tadd,sizeof(mtx_tadd));
    sysWrite8(fd, mtx_tffirst,sizeof(mtx_tffirst));
    sysWrite8(fd, mtx_textract,sizeof(mtx_textract));
    sysWrite8(fd, mtx_taddinv,sizeof(mtx_taddinv));
    sysWrite8(fd, mtx_tmultinv,sizeof(mtx_tmultinv));
    sysWrite8(fd, mtx_tnull,sizeof(mtx_tnull));
    sysWrite8(fd, mtx_tinsert,sizeof(mtx_tinsert));
    sysWrite32(fd,subfieldOrder,MTX_MAXSUBFIELDS);
    sysWrite8(fd,mtx_embed, MTX_MAXSUBFIELDS * MTX_MAXSUBFIELDORD);
    sysWrite8(fd,mtx_restrict, MTX_MAXSUBFIELDS * 256);

    fclose(fd);
    return(0);
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
