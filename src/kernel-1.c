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


int mtx_subfields[17];                  // public list of subfields, terminated with 0

static uint16_t minusone;                   // -1
static uint16_t *inc = NULL;                // a+1 = inc[a]
static uint16_t *FfFromIntTable = NULL;
static uint16_t *FfToIntTable = NULL;
static uint16_t subfieldsTable[17];     // internal list of subfields, terminated with 0xFFFF
static FEL *embeddingTables = NULL;         // combined embed/restrict tables

static unsigned short P = 0;         // Characteristic
static unsigned short Q = 0;         // Field order
static unsigned short Q1 = 0;        // Q-1, order of the multiplicative group
static unsigned short N;             // Degree over prime field, Q=P^N
static unsigned short Gen;           // Generator of the multiplicative group

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

static int loadEmbedAndRestrictTables(FILE *fd)
{
   uint16_t numberOfSubfields;
   sysRead16(fd, &numberOfSubfields, 1);
   MTX_ASSERT(numberOfSubfields <= 16);
        
   sysRead16(fd, subfieldsTable, numberOfSubfields);
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

   // Copy subfields to public table
   subfieldsTable[numberOfSubfields] = FF_INVALID;
   memset(mtx_subfields, 0, sizeof(mtx_subfields));
   for (int i = 0; i < numberOfSubfields && subfieldsTable[i] != FF_INVALID; ++i) {
      mtx_subfields[i] = (int) subfieldsTable[i];
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Tries loading the tables from "pXXXXX.zzz".
// Returns 1 on success or 0 if the file does not exist. If the file exists but cannot be opened
// or contains invalid data, the function raises an error, see @ref MtxError.
static int LoadTables_(int fieldOrder, const char* fileName)
{
   FILE* fd = sysFopen(fileName,"rb::lib:noerror");
   if (fd == NULL)
      return 0;

   // read header
   unsigned short info[5];
   if (fread((char *)info,sizeof(short),5,fd) != 5) {
      mtxAbort(MTX_HERE,"CANNOT READ TABLE HEADER");
   }

   P = info[1];
   Q = info[2];
   N = info[3];
   Gen = info[4];
   Q1 = Q - 1;
   ffOrder = (long) Q;
   ffChar = (long) P;
   ffGen = (Q == 2) ? 0 : 1;

   if (info[0] != MTX_ZZZVERSION) {
      mtxAbort(MTX_HERE,"Invalid table file: wrong version %d (expected %d)", info[0], MTX_ZZZVERSION);
   }

   if ((Q != fieldOrder) || (Q < 2) || (P < 2)
         || (P > Q) || (Q % P != 0)) {
      mtxAbort(MTX_HERE,"ERROR IN TABLE FILE HEADER");
   }


   // Read tables
   sysRead16(fd, &minusone, 1);
   inc = NREALLOC(inc, FEL, Q - 1);
   sysRead16(fd, inc, Q - 1);
   FfToIntTable = NREALLOC(FfToIntTable, uint16_t, Q);
   sysRead16(fd, FfToIntTable,Q);
   FfFromIntTable = NREALLOC(FfFromIntTable, uint16_t, Q);
   sysRead16(fd, FfFromIntTable,Q);
   loadEmbedAndRestrictTables(fd);

   fclose(fd);
   return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int LoadTables(int fieldOrder)
{
   static char fileName[50];
   snprintf(fileName, sizeof(fileName), "p%5.5d.zzz", fieldOrder);
   const int context = mtxBegin("Loading arithmetic tables: %s", fileName);
   int rc = LoadTables_(fieldOrder, fileName);
   mtxEnd(context);
   return rc;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Sets the field order.
/// This function sets the current field to GF(@em field) and initializes the field arithmetic.
/// Most kernel functions require that a field has been selected before they are used.
/// @param field Field order.
/// @return 0 on success, -1 otherwise.

int ffSetField(int field)
{
   MTX_ASSERT(sizeof(FEL) == 2);
   MTX_ASSERT(sizeof(unsigned short) == 2);
   MTX_ASSERT(sizeof(unsigned int) >= 4);

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
   MTX_ASSERT(noc >= 0);
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

size_t ffRowSizeUsed(int noc)
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

int ffFindPivot(PTR row,FEL *mark, int noc)
{
   register long i;
   register PTR p = row;

   for (i = 0; i < noc; ++i, ++p) {
      if (*p != FF_ZERO) {
         *mark = *p;
         return i;
      }
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

PTR ffAddRow(PTR dest, PTR src, int noc)
{
   register int i;
   register FEL *p1 = dest;
   register FEL *p2 = src;

   for (i = noc; i != 0; --i) {
      *p1 = ffAdd(*p1,*p2++);
      ++p1;
   }
   return dest;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void ffAddRowPartial(PTR dest, PTR src, int first, int noc)
{
   MTX_ASSERT(first >= 0 && first < noc);
   long i;

   PTR p1 = dest + first;
   PTR p2 = src + first;
   for (i = noc - first; i != 0; --i) {
      *p1 = ffAdd(*p1,*p2);
      p1++;
      p2++;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void ffMulRow(PTR row, FEL mark, int noc)
{
   register FEL *m;
   register long i;

   CHECKFEL(mark);
   if (mark == FF_ZERO) {
      m = row;
      for (i = noc; i != 0; --i) {
         *m++ = FF_ZERO;
      }

      // Fill unused space with zeroes.
      int rem = (noc * sizeof(FEL)) % sizeof(long);
      if (rem > 0)
         memset(m, 0, sizeof(long) - rem);
   } else {
      m = row;
      for (i = noc; i != 0; --i) {
         if (*m != FF_ZERO) {
            unsigned int x = *m + mark;
            *m = (x >= Q1) ? x - Q1 : x;
         }
         ++m;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void ffAddMulRow(PTR row1, PTR row2, FEL f, int noc)
{
   CHECKFEL(f);
   if (f == FF_ONE) {
      ffAddRow(row1, row2, noc);
      return;
   }
   else if (f != FF_ZERO) {
      FEL* p1 = row1;
      FEL* p2 = row2;
      for (int i = noc; i > 0; --i) {
         *p1 = ffAdd(*p1,ffMul(*p2,f));
         ++p1;
         ++p2;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void ffAddMulRowPartial(PTR dest, PTR src, FEL f, int firstcol, int noc)
{
   CHECKFEL(f);
   MTX_ASSERT(firstcol >= 0 && firstcol < noc);

   if (f == FF_ONE) {
      ffAddRowPartial(dest, src, firstcol, noc);
      return;
   }
   else if (f != FF_ZERO) {
      PTR p1 = dest + firstcol;
      PTR p2 = src + firstcol;
      for (int i = noc - firstcol; i != 0; --i) {
         if (*p2 != FF_ZERO) {
            *p1 = ffAdd(*p1,ffMul(*p2,f));
         }
         ++p1;
         ++p2;
      }
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

void ffMapRow(PTR row, PTR matrix, int nor, int noc, PTR result)
{
   ffMulRow(result, FF_ZERO, noc);

   FEL *brow = row;
   FEL *m = matrix;
   for (int i = nor; i > 0; --i) {
      FEL f = *brow++;

      if (f != FF_ZERO) {
         FEL *v = m;
         FEL *r = result;
         if (f == FF_ONE) {
            for (int k = noc; k != 0; --k) {
               *r = ffAdd(*r,*v);
               ++r;
               ++v;
            }
         } else {
            for (int k = noc; k != 0; --k) {
               *r = ffAdd(*r,ffMul(*v,f));
               ++r;
               ++v;
            }
         }
      }

      m += ffRowSize(noc) / sizeof(FEL);   /* next row */
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

FEL ffScalarProduct(PTR a, PTR b, int noc)
{
   FEL *ap = a;
   FEL *bp = b;
   FEL f = FF_ZERO;
   for (int i = noc; i > 0; --i) {
      f = ffAdd(f,ffMul(*ap++,*bp++));
   }
   return f;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void ffExtractColumn(PTR mat, int nor, int noc, int col, PTR result)
{
   FEL *x = mat + col;
   FEL *y = result;

   const size_t step = ffRowSize(noc) / sizeof(FEL);
   for (int count = nor; count > 0; --count) {
      *y++ = *x;
      x += step;
   }
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
