////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Finite field arithmetic and common functions
// This is the "Large" version for field orders q <= 65535.
//
// (C) Copyright 1998-2014 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "maketab-1.h"
#include <stdlib.h>
#include <string.h>



/* ------------------------------------------------------------------
   Data
   ------------------------------------------------------------------ */

MTX_DEFINE_FILE_INFO

FEL minusone;			/* -1 */
FEL *inc = NULL;		/* a+1 = inc[a] */

static unsigned short P = 0;		/* Characteristic */
static unsigned short Q = 0;		/* Field order */
static unsigned short Q1 = 0;		/* Q-1 */
static unsigned short N;
static unsigned short Gen;
static unsigned short ppwr[MAXPWR];
static unsigned short ppindex[MAXPWR];
static unsigned short poly[MAXPWR];

static long IPR = 0;            /* No. of long ints per row */

/* ------------------------------------------------------------------
   Argument checking macros
   ------------------------------------------------------------------ */

#if defined(DEBUG)

#define CHECKRANGE(x,lo,hi) if ((x)<(lo)||(x)>(hi)) {\
	fprintf(stderr,"%ld <= %ld <= %ld ?\n",(long)(lo),\
		(long)(x),(long)(hi));\
	MTX_ERROR("RANGE CHECK ERROR");}
#define CHECKFILE(x)  CHECKRANGE(file,0,MAXFILES)
#define CHECKCOL(x)  CHECKRANGE(x,1,FfNoc)
#define CHECKFEL(x) { \
	if ((x) != 0xFFFF && ((x) > Q-2)) \
		MTX_ERROR("range check error"); \
	}

#else

#define CHECKRANGE(x,lo,hi)
#define CHECKCOL(x)
#define CHECKFILE(x)
#define CHECKFEL(x)

#endif


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Set the row length.
/// This function sets the current row size, which is used for low-level row operations
/// such as FfAddRow().
/// @param noc Number of columns.
/// @return 0 on success, -1 otherwise.

int FfSetNoc(int ncols)
{
   FfNoc = ncols;
   FfCurrentRowSizeIo = FfNoc * sizeof(FEL);

   IPR = (FfNoc*sizeof(FEL) + (long) - 11) / (sizeof(long));
   FfCurrentRowSize = IPR * sizeof(long);
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Sets the field order.
/// This function sets the current field to GF(@em field) and initializes the field arithmetic.
/// Most kernel functions require that a field has been selected before they are used.
/// @param field Field order.
/// @return 0 on success, -1 otherwise.

int FfSetField(int field)
{
   static char filename[50];
   FILE *fd;
   unsigned short info[5];

   if (field != FfOrder)
   {
       Ff1MakeTables(field);

      sprintf(filename,"p%5.5d.zzz",field);
      if ((fd = SysFopen(filename,FM_READ|FM_LIB)) == NULL ||
	    fread((char *)info,sizeof(short),5,fd) != 5)
      {	perror(filename);
	 MTX_ERROR("Error opening table file");
	 return -1;
      }
      P = info[1];
      Q = info[2];
      N = info[3];
      Q1 = Q - 1;
      Gen = info[4];
      FfOrder = (long) Q;
      FfChar = (long) P;
      if (Q == 2) FfGen = 0; else FfGen = 1;
      if ((info[0]) != ZZZVERSION || N >= MAXPWR || Q < 2 || P < 2 || P > Q || Q % P != 0) {
	 MTX_ERROR("ERROR IN TABLE FILE HEADER");
	 return -1;
      }
      if (inc != NULL) free(inc);
      inc = NALLOC(FEL,Q-1);
      if (
	    fread(poly,sizeof(short),(size_t) N+1,fd) != (size_t)N+1 ||
	    fread(ppwr,sizeof(short),(size_t) N,fd) != (size_t)N ||
	    fread(ppindex,sizeof(short),(size_t) N,fd) != (size_t)N ||
	    fread(&minusone,sizeof(short),1,fd) != 1 ||
	    fread(inc,sizeof(short),(size_t)(Q-1),fd) != (int) Q-1
	 )
      {	perror(filename);
	 exit(EXIT_ERR);
      }
      if (poly[N] != 1) MTX_ERROR("ERROR IN TABLE FILE");
      fclose(fd);
   }
   return 0;
}


/* ------------------------------------------------------------------
   FfAdd(), FfSub(), FfMul(), FfDiv() - Field arithmetic
   ------------------------------------------------------------------ */

FEL FfAdd(FEL a, FEL b)
{	register FEL x;

   CHECKFEL(a);
   CHECKFEL(b);

   if (b == FF_ZERO) return a;
   if (a == FF_ZERO) return b;
   if (a >= b)
   {	x = inc[a-b];
      if (x == FF_ZERO) return FF_ZERO;
      if ((x += b) >= Q1) x -= Q1;
      return x;
   }
   x = inc[b-a];
   if (x == FF_ZERO) return FF_ZERO;
   if ((x += a) >= Q1) x -= Q1;
   return x;
}


FEL FfSub(FEL a, FEL b)
{	register FEL x, bb;

   CHECKFEL(a);
   CHECKFEL(b);

   if (b == FF_ZERO) return a;
   if (b == a) return FF_ZERO;
   if ((bb = b + minusone) >= Q1) bb -= Q1;	/* bb = -b */
   if (a == FF_ZERO) return bb;
   if (a >= bb)
   {	x = inc[a-bb];
      if ((x += bb) >= Q1) x -= Q1;
      return x;
   }
   x = inc[bb-a];
   if ((x += a) >= Q1) x -= Q1;
   return x;
}


FEL FfMul(FEL a, FEL b)
{	register FEL c;

   CHECKFEL(a);
   CHECKFEL(b);

   if (b == FF_ZERO || a == FF_ZERO) return FF_ZERO;
   if ((c = a + b) >= Q1) c -= Q1;
   return (c);
}


FEL FfDiv(FEL a, FEL b)
{
   CHECKFEL(a);
   CHECKFEL(b);
   if (b == FF_ZERO) {
      MTX_ERROR("Division by zero");
      return FF_ZERO;
   }
   if (a == FF_ZERO) return FF_ZERO;
   if (a >= b)
      return (a-b);
   else
      return ((a+Q1)-b);
}

FEL FfNeg(FEL a)
{
   return FfSub(FF_ZERO, a);
}


FEL FfInv(FEL a)
{
   return FfDiv(FF_ONE, a);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Calculate row size.
/// Returns the number of bytes occupied in memory by a row of @em noc Elements.
/// The row size is always a multiple of <tt>sizeof(long)</tt>. Depending on the number of
/// columns there may be unused padding bytes at the end of the row.

size_t FfRowSize(int noc)
{
   if (noc == 0) {
      return 0;
   }
   return ((noc * sizeof(FEL) + sizeof(long) - 1) / sizeof(long)) * sizeof(long);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Number of used bytes in a row.
/// This function returns the number of bytes that are actually used by a row of @em noc Elements,
/// i.e., without counting the padding bytes. This number is less than or equal to
/// <tt>FfRowSize(noc)</tt>.

size_t FfTrueRowSize(int noc)
{
   return noc * sizeof(FEL);
}


/* ------------------------------------------------------------------
   FfEmbed() - Embed a subfield element.

   return value: Embedded field element.
   ------------------------------------------------------------------ */

FEL FfRestrict(FEL a, int subfield)
{
   int i;
   long l;

   /* If a is zero, return zero.
      -------------------------- */
   if (a == FF_ZERO) return FF_ZERO;

   /* If a is non-zero, find out the degree of the subfield
      relative to the current field, and divide by this
      number.
      ------------------------------------------------------ */
   for (l = subfield, i = 1; l < Q; l *= subfield, ++i);
   if (l != Q)
      MTX_ERROR("ILLEGAL SUBFIELD");
   return (a / i);
}


/* ------------------------------------------------------------------
   FfEmbed() - Embed a subfield element.

   return value: Embedded field element.
   ------------------------------------------------------------------ */

FEL FfEmbed(FEL a, int subfield)

{
   int i;
   long l;

   /* If a is zero, return zero.
      -------------------------- */
   if (a == FF_ZERO) return FF_ZERO;

   /* If a is non-zero, find out the degree of the subfield
      relative to the current field, and multiply by this
      number.
      ------------------------------------------------------ */
   for (l = subfield, i = 1; l < Q; l *= subfield, ++i);
   if (l != Q)
      MTX_ERROR("ILLEGAL SUBFIELD");
   return (a * i);
}





/* ------------------------------------------------------------------
   FfFindPivot()  -  Find first non-zero mark in row

   return value: column (0 if row is zero)
   ------------------------------------------------------------------ */

int FfFindPivot(PTR row,FEL *mark)

{
   register long i;
   register PTR p = row;

   for (i = 1; i <= FfNoc; ++i, ++p)
   {	if (*p != FF_ZERO)
      {	*mark = *p;
	 return (i);
      }
   }
   return 0;
}




////////////////////////////////////////////////////////////////////////////////////////////////////
/// Add two rows.
/// This function adds src to dest. Field order and row size must have been set before.
/// @param dest The row to add to.
/// @param src The row to add.
/// @return Always returns dest.

PTR FfAddRow(PTR dest, PTR src)
{
   register int i;
   register FEL *p1 = dest;
   register FEL *p2 = src;

   for (i = FfNoc; i != 0; --i)
   {
      *p1 = FfAdd(*p1,*p2++);
      ++p1;
   }
   return dest;
}


/* ------------------------------------------------------------------
   FfMulRow()  -  multiply row by field element

   return value: none
   ------------------------------------------------------------------ */

void FfMulRow(PTR row, FEL mark)

{	register FEL *m;
	register long i;

	CHECKFEL(mark);
	if (mark == FF_ZERO)
	{	m = row;
		for (i = FfNoc; i != 0; --i)
		{	*m++ = FF_ZERO;
		}
	}
	else
	{	m = row;
		for (i = FfNoc; i != 0; --i)
		{	if (*m != FF_ZERO)
			{	if ((*m += mark) >= Q1)
					*m -= Q1;
			}
			++m;
		}
	}
}


/* ------------------------------------------------------------------
   FfAddMulRow() - row1 += f * row2

   return value: none
   ------------------------------------------------------------------ */

void FfAddMulRow(PTR row1, PTR row2, FEL f)

{	register long i;
	register FEL *p1, *p2;

	CHECKFEL(f);
	if (f == FF_ZERO) return;
	if (f == FF_ONE)
	{	FfAddRow(row1,row2);
		return;
	}
	p1 = row1;
	p2 = row2;
	for (i = FfNoc; i != 0; --i)
	{	*p1 = FfAdd(*p1,FfMul(*p2,f));
		++p1;
		++p2;
	}
}



/* ------------------------------------------------------------------
   FfFromInt(), FfToInt() - Convert between long int and FEL
   ------------------------------------------------------------------ */

static void polmul(a, b)
unsigned short *a, *b;

{	unsigned short x[2*MAXPWR], f;
	int i, k;

	memset(x,0,sizeof(x));
	for (i = 0; i < (int) N; ++i)
		for (k = 0; k < (int)N; ++k)
			x[i+k] = (unsigned short)
				((x[i+k] + (long)a[i]*b[k]) % P);
	for (i = 2*(int)N-2; i >= (int)N; --i)
	{	if ((f = x[i]) == 0) continue;
		f = (unsigned short) P - f;
		for (k = (int)N; k >= 0; --k)
		{	x[i-(int)N+k] = (unsigned short)
			    ((x[i-(int)N+k] + (long) f * poly[k]) % P);
		}
	}
	memcpy(a,x,(size_t)N * sizeof(short));
}






FEL FfFromInt(int l)

{	register FEL f = FF_ZERO;
	register int i;

	l = l % Q;
	if (l < 0) l += Q;
	for (i = (int)N-1; i >= 0; --i)
	{	while (l >= ppwr[i])
		{	f = FfAdd(f,ppindex[i]);
			l -= ppwr[i];
		}
	}
	return f;
}


int FfToInt(FEL f)

{	unsigned short i;
	int l, m;
	unsigned short a[MAXPWR], b[MAXPWR];

	if (f == FF_ZERO) return 0;
	if (N == 1)
	{	m = Gen;
		l = 1;
		for (i = 1; f != 0; i <<= 1)
		{	if ((f & i) != 0)
			{	l = (l * m) % P;
				f &= ~i;
			}
			m = (m * m) % P;
		}
	}
	else
	{	memset(a,0,sizeof(a));
		a[0] = 1;			/* a(x) = 1 */
		memset(b,0,sizeof(b));
		b[1] = 1;			/* b(x) = x */
		for (i = 1; f != 0; i <<= 1)
		{	if ((f & i) != 0)
			{	polmul(a,b);
				f &= ~i;
			}
			polmul(b,b);
		}
		l = 0;
		m = 1;
		for (i = 0; i < (short) N; ++i)
		{	l += a[i] * m;
			m *= P;
		}
	}
	return l;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Multiply a vector by a matrix.
/// This function multiplies the vector @em row from the right by the matrix @em mat and
/// stores the result into @em result.
/// The number of columns in both @em mat and @em result is determined by the current row size.
/// (see FfNoc()).
/// @attention @em result and @em row must not overlap. Otherwise the result is
/// undefined.
/// @param row The source vector (FfNoc columns).
/// @param matrix The matrix (nor by FfNoc).
/// @param nor number of rows in the matrix.
/// @param[out] result The resulting vector (@em nor columns).

void FfMapRow(PTR row, PTR matrix, int nor, PTR result)
{
   memset(result, 0, FfCurrentRowSize);

   FEL *brow = row;
   for(int i = nor; i > 0; --i) {
       FEL *m = matrix;
       FEL f = *brow++;

       if (f != FF_ZERO) {
	   FEL *v = m;
	   FEL *r = result;
	   int k = FfNoc;
	   if (f == FF_ONE) {
	       for (; k != 0; --k) {
		   *r = FfAdd(*r,*v);
		   ++r;
	       }
	   } else {
	       for (; k != 0; --k) {
		   *r = FfAdd(*r,FfMul(*v,*r));
		   ++r;
		   ++v;
	       }
	   }
       }

       m += FfCurrentRowSize;                 /* next row */
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Scalar Product of Two Vectors.
/// Given two vectors @f$a=(a_i)@f$ and @f$b=(b_i)@f$, this function calculates the
/// scalar product @f$p=\sum_ia_ib_i@f$.
/// @param a The first vector.
/// @param b The second vector.
/// @return Scalar product of the two vectors.

FEL FfScalarProduct(PTR a, PTR b)
{
   FEL *ap = a;
   FEL *bp = b;
   FEL f = FF_ZERO;
   for (int i = FfNoc; i >= 0; --i) {
      f = FfAdd(f,FfMul(*ap++,*bp++));
   }
   return f;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Extract one column of a matrix.
/// This function extracts one column out of a matrix and converts it into a row vector.
/// The number of columns of the matrix must have been set with FfSetNoc(), and the number of rows
/// must be passed in the call. The result is a row with @a nor entries, and the output buffer,
/// @a result, must be large enough to store a vector of this size.
/// @a mat and @a result must not overlap. If they do, the result is undefined.
///
/// @param mat Pointer to the matrix.
/// @param nor Number of rows in matrix.
/// @param col Column to extract (starting with 0).
/// @param result Pointer to buffer for the extracted column.

void FfExtractColumn(PTR mat, int nor, int col, PTR result)
{
   FEL *x = mat + col;
   FEL *y = result;

   for (int count = nor; count > 0; --count) {
      *y++ = *x;
      x += FfCurrentRowSize;
   }
}

