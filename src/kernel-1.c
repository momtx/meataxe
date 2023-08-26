////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Finite field arithmetic and common functions
// This is the "Large" version for field orders q <= 65535.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>
#include <string.h>

// The "large" arithmetics kernel supports fiels order up to 63001.
//
// The internal data representation is based on a generator for the multiplicative group which is
// chosen for each field as follows. If q = p is prime, the generator is the smallest natural
// number g which generates the multiplicative group {1, 2, 3, 4, ... p-1} in ℤₚ.
// If q = p^n is not prime, we construct GF(q) = F[x]/(p(x)), where p(x) is the Conway
// polynomial associated with q (see the table in maketab-1.c). The generator, in this case,
// is g = x + (p(x)).
// Nonzero field elements are internally represented as their logarithm base g, i.e., the
// unit element is represented as 0, the generator g as 1, and so on. The zero element is
// represented by the special value 65536, which can never occur as representatzion of a
// nonzero field element.
// Note that the prime field elements are NOT represented by 0..p-1, even for prime fields.
// All calculations must be performed with FfAdd/FfMul. Always use the constants
// FF_ZERO (=0xFFFF) and FF_ONE (=0) for the zero and unit element.
//
// External representation (defined by @ref FfToInt and @ref FfFromInt):
// For prime fields, the external representation are the integers {0,1, ..., p-1} that were
// used for the field construction, and standard integer arithmetic modulo p can be used for
// calculations.
// For non-prime fields, the field elements are represented by the integer numbers {0..,q-1},
// which are found by inserting p into the element's representative polynomial. The numbers
// {0..p} correspond to the prime field.



FEL minusone;                   /* -1 */
FEL *inc = NULL;                /* a+1 = inc[a] */
FEL *FfFromIntTable = NULL;
unsigned short *FfToIntTable = NULL;
FEL subfieldsTable[17];         // list of proper subfields, terminated with 0xFFFF
FEL *embeddingTables = NULL;    // combined embed/restrict tables

static unsigned short P = 0;            // Characteristic
static unsigned short Q = 0;            // Field order
static unsigned short Q1 = 0;           // Q-1, order of the multiplicative group
static unsigned short N;                // Degree over prime field, Q=P^N
static unsigned short Gen;              // Generator of the multiplicative group

static long IPR = 0;            /* No. of long ints per row */

#define FF_INVALID 0xFFFE

////////////////////////////////////////////////////////////////////////////////////////////////////
// Argument checking macros
////////////////////////////////////////////////////////////////////////////////////////////////////

#define ISFELX(x, q) ((x) == FF_ZERO || (x) < (q) - 1)

#if defined(MTX_DEBUG)

#define CHECKRANGE(x,lo,hi) if ((x) < (lo) || (x) > (hi)) { \
      fprintf(stderr,"%ld <= %ld <= %ld ?\n",(long)(lo), \
              (long)(x),(long)(hi)); \
      mtxAbort(MTX_HERE,"RANGE CHECK ERROR");}
#define CHECKFILE(x)  CHECKRANGE(file,0,MAXFILES)
#define CHECKCOL(x)  CHECKRANGE(x,1,ffNoc)
#define CHECKFEL(x) { \
      if ((x) != 0xFFFF && ((x) > Q - 2)) \
      mtxAbort(MTX_HERE,"range check error"); \
}

#else

#define CHECKRANGE(x,lo,hi)
#define CHECKCOL(x)
#define CHECKFILE(x)
#define CHECKFEL(x)

#endif

////////////////////////////////////////////////////////////////////////////////////////////////////

void ffSetNoc(int ncols)
{
   ffNoc = ncols;
   IPR = (ffNoc * sizeof(FEL) + sizeof(long) - 1) / (sizeof(long));
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int loadEmbedAndRestrictTables(FILE *fd)
{
   FEL numberOfSubfields;
   if (fread(&numberOfSubfields, sizeof(FEL), 1, fd) != 1 || numberOfSubfields > 16) {
      mtxAbort(MTX_HERE,"Currupt table file (number of subfields)");
      return -1;
   }
   if (fread(subfieldsTable, sizeof(FEL), numberOfSubfields, fd) != numberOfSubfields) {
      mtxAbort(MTX_HERE,"Currupt table file (subfield orders)");
      return -1;
   }
   subfieldsTable[numberOfSubfields] = FF_INVALID;

   size_t tblSize = 0;
   for (FEL* sf = subfieldsTable; *sf != FF_INVALID; ++sf) {
      if (*sf >= Q) {
         mtxAbort(MTX_HERE,"Corrupt table file (subfield order)");
      }
      tblSize += *sf + Q;
   }
   embeddingTables = NREALLOC(embeddingTables, FEL, tblSize);
   if (fread(embeddingTables, sizeof(FEL), tblSize, fd) != tblSize) {
      mtxAbort(MTX_HERE,"Corrupt table file (subfield embeddings)");
      return -1;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Tries loading the tables from "pXXXXX.zzz".
// Returns 1 on success or 0 if the file does not exist. If the file exists but cannot be opened
// or contains invalid data, the function raises an error, see @ref MtxError.
static int LoadTables(int fieldOrder)
{
   static char filename[50];
   sprintf(filename,"p%5.5d.zzz",fieldOrder);
   FILE* fd = sysFopen(filename,"rb::lib:noerror");
   int ok = fd != NULL;

   // read header
   unsigned short info[5];
   if (ok && (fread((char *)info,sizeof(short),5,fd) != 5)) {
      mtxAbort(MTX_HERE,"CANNOT READ TABLE HEADER");
      ok = 0;
   }
   if (ok) {
      P = info[1];
      Q = info[2];
      N = info[3];
      Gen = info[4];
      Q1 = Q - 1;
      ffOrder = (long) Q;
      ffChar = (long) P;
      ffGen = (Q == 2) ? 0 : 1;

      if (info[0] != ZZZVERSION) {
         mtxAbort(MTX_HERE,"Invalid table file: wrong version %d (expected %d)", info[0], ZZZVERSION);
         ok = 0;
      }

      if ((Q != fieldOrder) || (Q < 2) || (P < 2)
          || (P > Q) || (Q % P != 0)) {
         mtxAbort(MTX_HERE,"ERROR IN TABLE FILE HEADER");
         ok = 0;
      }
   }

   // Read tables
   if (ok) {
      inc = NREALLOC(inc, FEL, Q - 1);
      FfFromIntTable = NREALLOC(FfFromIntTable, FEL, Q);
      FfToIntTable = NREALLOC(FfToIntTable, unsigned short, Q);
      if ((fread(&minusone,sizeof(short),1,fd) != 1) ||
          (fread(inc,sizeof(short),(size_t)(Q - 1),fd) != (size_t) Q - 1) ||
          (fread(FfToIntTable,sizeof(short),(size_t)Q,fd) != (size_t) Q) ||
          (fread(FfFromIntTable,sizeof(short),(size_t)Q,fd) != (size_t) Q) ||
          loadEmbedAndRestrictTables(fd) != 0 ||
          ferror(fd)
          ) {
         mtxAbort(MTX_HERE,"COULD NOT LOAD ARITHMETIC TABLE FILE");
         ok = 0;
      }
   }

   if (fd != NULL) {
      fclose(fd);
   }
   return ok;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Sets the field order.
/// This function sets the current field to GF(@em field) and initializes the field arithmetic.
/// Most kernel functions require that a field has been selected before they are used.
/// @param field Field order.
/// @return 0 on success, -1 otherwise.

int ffSetField(int field)
{
   MTX_ASSERT(sizeof(FEL) == 2, -1);
   MTX_ASSERT(sizeof(unsigned short) == 2, -1);
   MTX_ASSERT(sizeof(unsigned int) >= 4, -1);

   if ((field == ffOrder) || (field < 2)) {
      return 0;
   }
   if (!LoadTables(field)) {
      ffMakeTables(field);
      if (!LoadTables(field)) {
         mtxAbort(MTX_HERE,"COULD NOT LOAD ARITHMETIC TABLE FILE");
         return -1;
      }
   }

   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

FEL ffAdd(FEL a, FEL b)
{
   CHECKFEL(a);
   CHECKFEL(b);

   if (b == FF_ZERO) {return a;}
   if (a == FF_ZERO) {return b;}
   unsigned int x;
   if (a >= b) {
      // a + b = (a/b + 1) * b
      x = inc[a - b];
      if (x == FF_ZERO) {return FF_ZERO;}
      if ((x += b) >= Q1) {x -= Q1;}
   } else {
      // a + b = (b/a + 1) * a
      x = inc[b - a];
      if (x == FF_ZERO) {return FF_ZERO;}
      if ((x += a) >= Q1) {x -= Q1;}
   }
   return x;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

FEL ffSub(FEL a, FEL b)
{
   CHECKFEL(a);
   CHECKFEL(b);

   if (b == FF_ZERO) {return a;}
   if (b == a) {return FF_ZERO;}

   unsigned int minusb = (unsigned int) b + minusone;
   if (minusb >= Q1) {minusb -= Q1;}
   if (a == FF_ZERO) {return minusb;}

   // Same as ffAdd(a,-b), but we can omit the check for FF_ZERO since a != b.
   unsigned int x;
   if (a >= minusb) {
      // a - b = (a/(-b) + 1) * (-b);
      x = inc[a - minusb] + minusb;
   } else {
      // a - b = ((-b)/a + 1) * a;
      x = inc[minusb - a] + a;
   }
   return x >= Q1 ? (x - Q1) : x;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

FEL ffMul(FEL a, FEL b)
{
   CHECKFEL(a);
   CHECKFEL(b);

   if ((b == FF_ZERO) || (a == FF_ZERO)) {return FF_ZERO;}
   unsigned c = a + b;
   return c >= Q1 ? (c - Q1) : c;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

FEL ffDiv(FEL a, FEL b)
{
   CHECKFEL(a);
   CHECKFEL(b);
   if (b == FF_ZERO) {
      mtxAbort(MTX_HERE,"Division by zero");
      return FF_ZERO;
   }
   if (a == FF_ZERO) {return FF_ZERO;}
   if (a >= b) {
      return a - b;
   } else {
      return Q1 - (b - a);
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

FEL ffNeg(FEL a)
{
   return ffSub(FF_ZERO, a);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

FEL ffInv(FEL a)
{
   CHECKFEL(a);
   if (a == FF_ZERO) {
      mtxAbort(MTX_HERE,"Division by zero");
      return FF_ZERO;
   }
   if (a == FF_ONE) {
      return FF_ONE;
   }
   return Q1 - a;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

size_t ffRowSize(int noc)
{
   MTX_ASSERT(noc >= 0, 0);
   if (noc == 0) {
      return 0;
   }
   return ((noc * sizeof(FEL) + sizeof(long) - 1) / sizeof(long)) * sizeof(long);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

ssize_t ffSize(int nor, int noc)
{
   return nor == 0 ? 0 : nor * (ssize_t) ffRowSize(noc);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Number of used bytes in a row.
/// This function returns the number of bytes that are actually used by a row of @em noc Elements,
/// i.e., without counting the padding bytes. This number is less than or equal to
/// <tt>ffRowSize(noc)</tt>.

size_t ffTrueRowSize(int noc)
{
   return noc * sizeof(FEL);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Embed/Restrict tables: for each subfield
//  - FEL embed[s] - embedding into F(q)
//  - FEL restrict[q] - restriction to F(s) or 0xFFFE 
// all tables are contiguous without padding. The zero element is handled in code and not
// included in the tables
// Subfield table: 
//  - uint16_t[]: list of subfield order, terminated with 0xFFFF
//    e.g., for q=3^6:   [3, 9, 27, 0xFFFF]


// Returns a pointer to the combined embed/restrict table for F(r)<F(q).
// The table has size r + q.
static FEL* getEmbeddingTable(int r)
 {
   // Look up the subfield for every call, assuming ffEmbed() performance is not critical.
   FEL* sptr = subfieldsTable;
   FEL* tptr = embeddingTables;
   while (*sptr != r) {
      if (*sptr == 0xFFFF) {
         mtxAbort(MTX_HERE,"Bad subfield. Cannot embed F(%d) into F(%d).", r, Q);
         return NULL;
      }
      if (*sptr >= Q) {
         mtxAbort(MTX_HERE,"Corrupt subfield table.");
         return NULL;
      }
      tptr += *sptr + Q;
      ++sptr;
   }
   return tptr;
}
 

FEL ffRestrict(FEL a, int subfield)
{
   if (a == FF_ZERO) {return FF_ZERO;}
   if (a == FF_ONE) {return FF_ONE;}
   CHECKFEL(a);

   FEL* table = getEmbeddingTable(subfield);
   if (table == NULL)
      return FF_ZERO;
   FEL result = table[subfield + a];
   if (result == FF_INVALID) {
      mtxAbort(MTX_HERE,"%s(): Element %u is not in subfield F(%u).", __func__, a, subfield);
      return FF_ZERO;
   }
   return result;
}
 
////////////////////////////////////////////////////////////////////////////////////////////////////

FEL ffEmbed(FEL a, int subfield)
{
   if (a == FF_ZERO) {return FF_ZERO;}
   if (a == FF_ONE) {return FF_ONE;}

   if (!ISFELX(a, subfield)) {
      mtxAbort(MTX_HERE,"FfEmbed: subfield element 0x%x not in F(%u)", a, subfield);
      return FF_ZERO;
   }
   FEL* table = getEmbeddingTable(subfield);
   return table ? table[a] : FF_ZERO;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

int ffFindPivot(PTR row,FEL *mark)
{
   register long i;
   register PTR p = row;

   for (i = 0; i < ffNoc; ++i, ++p) {
      if (*p != FF_ZERO) {
         *mark = *p;
         return i;
      }
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

PTR ffAddRow(PTR dest, PTR src)
{
   register int i;
   register FEL *p1 = dest;
   register FEL *p2 = src;

   for (i = ffNoc; i != 0; --i) {
      *p1 = ffAdd(*p1,*p2++);
      ++p1;
   }
   return dest;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

PTR ffAddRowPartial(PTR dest, PTR src, int first)
{
   MTX_ASSERT(first >= 0 && first < ffNoc, dest);
   long i;

   PTR p1 = dest + first;
   PTR p2 = src + first;
   for (i = ffNoc - first; i != 0; --i) {
      *p1 = ffAdd(*p1,*p2);
      p1++;
      p2++;
   }
   return dest;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

void ffMulRow(PTR row, FEL mark)
{
   register FEL *m;
   register long i;

   CHECKFEL(mark);
   if (mark == FF_ZERO) {
      m = row;
      for (i = ffNoc; i != 0; --i) {
         *m++ = FF_ZERO;
      }

      // Fill unused space with zeroes.
      int rem = (ffNoc * sizeof(FEL)) % sizeof(long);
      if (rem > 0)
         memset(m, 0, sizeof(long) - rem);
   } else {
      m = row;
      for (i = ffNoc; i != 0; --i) {
         if (*m != FF_ZERO) {
            unsigned int x = *m + mark;
            *m = (x >= Q1) ? x - Q1 : x;
         }
         ++m;
      }
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

void ffAddMulRow(PTR row1, PTR row2, FEL f)
{
   register long i;
   register FEL *p1, *p2;

   CHECKFEL(f);
   if (f == FF_ZERO) {return;}
   if (f == FF_ONE) {
      ffAddRow(row1,row2);
      return;
   }
   p1 = row1;
   p2 = row2;
   for (i = ffNoc; i != 0; --i) {
      *p1 = ffAdd(*p1,ffMul(*p2,f));
      ++p1;
      ++p2;
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

void ffAddMulRowPartial(PTR dest, PTR src, FEL f, int firstcol)
{
   CHECKFEL(f);
   MTX_ASSERT(firstcol >= 0 && firstcol < ffNoc,);

   if (f == FF_ZERO) {
      return;
   }
   if (f == FF_ONE) {
      ffAddRowPartial(dest, src, firstcol);
      return;
   }
   PTR p1 = dest + firstcol;
   PTR p2 = src + firstcol;
   for (int i = ffNoc - firstcol; i != 0; --i) {
      if (*p2 != FF_ZERO) {
         *p1 = ffAdd(*p1,ffMul(*p2,f));
      }
      ++p1;
      ++p2;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

FEL ffFromInt(const int l)
{
   int m = l % Q;
   if (m < 0) { m += Q;}
   return FfFromIntTable[m];
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int ffToInt(FEL f)
{
   CHECKFEL(f);
   return (f == FF_ZERO) ? 0 : FfToIntTable[f];
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void ffMapRow(PTR row, PTR matrix, int nor, PTR result)
{
   ffMulRow(result, FF_ZERO);

   FEL *brow = row;
   FEL *m = matrix;
   for (int i = nor; i > 0; --i) {
      FEL f = *brow++;

      if (f != FF_ZERO) {
         FEL *v = m;
         FEL *r = result;
         int k = ffNoc;
         if (f == FF_ONE) {
            for (; k != 0; --k) {
               *r = ffAdd(*r,*v);
               ++r;
               ++v;
            }
         } else {
            for (; k != 0; --k) {
               *r = ffAdd(*r,ffMul(*v,f));
               ++r;
               ++v;
            }
         }
      }

      m += ffRowSize(ffNoc) / sizeof(FEL);   /* next row */
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Scalar Product of Two Vectors.
/// Given two vectors @f$a=(a_i)@f$ and @f$b=(b_i)@f$, this function calculates the
/// scalar product @f$p=\sum_ia_ib_i@f$.
/// @param a The first vector.
/// @param b The second vector.
/// @return Scalar product of the two vectors.

FEL ffScalarProduct(PTR a, PTR b)
{
   FEL *ap = a;
   FEL *bp = b;
   FEL f = FF_ZERO;
   for (int i = ffNoc; i > 0; --i) {
      f = ffAdd(f,ffMul(*ap++,*bp++));
   }
   return f;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Extract one column of a matrix.
/// This function extracts one column out of a matrix and converts it into a row vector.
/// The number of columns of the matrix must have been set with ffSetNoc(), and the number of rows
/// must be passed in the call. The result is a row with @a nor entries, and the output buffer,
/// @a result, must be large enough to store a vector of this size.
/// @a mat and @a result must not overlap. If they do, the result is undefined.
///
/// @param mat Pointer to the matrix.
/// @param nor Number of rows in matrix.
/// @param col Column to extract (starting with 0).
/// @param result Pointer to buffer for the extracted column.

void ffExtractColumn(PTR mat, int nor, int col, PTR result)
{
   FEL *x = mat + col;
   FEL *y = result;

   const size_t step = ffRowSize(ffNoc) / sizeof(FEL);
   for (int count = nor; count > 0; --count) {
      *y++ = *x;
      x += step;
   }
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
