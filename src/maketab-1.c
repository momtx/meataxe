////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Calculate arithmetic tables (large fields, q < 65536)
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXPWR 16

// Polynomial in Z[p]
typedef uint16_t POLY[MAXPWR + 1];


////////////////////////////////////////////////////////////////////////////////////////////////////
// Data to be written into the table file,
////////////////////////////////////////////////////////////////////////////////////////////////////

static uint16_t* inc = NULL;              // Increment table
static uint16_t Minusone;

// Mapping of internal to external representation for nonzero field elements. This table
// has Q elements, but only positions 0...Q-2 are actually used. The zero element is handled
// separately.
static uint16_t* FfToIntTable = NULL;

// Mapping of external to internal representation.
// - The zero element is at position 0, i.e., FfFromIntTable[0] is always 0xFFFF.
// - For nonzero field elements, FfFromIntTable[i] is the logarithm base Gen of the
//   i-th field element. In particular, the unit element is always at position 1,
//   i.e., FfFromIntTable[1] is always 0.
static uint16_t* FfFromIntTable = NULL;

static uint16_t numberOfSubfields;
static uint16_t subfieldOrders[20];
static uint16_t* embeddingTables = NULL;        // Joined embed/restrict tables for all subfields
size_t embeddingTablesSize = 0;                 // Size in uint16_t units
static uint32_t P;        // Characteristic of the field
static uint32_t Q;        // Order of the field
static uint32_t N;        // Q = P^N
static uint16_t Gen;      // Generator

static FILE *fd;                // Output file

// The defining polynomial d(x). For prime fields, the polynomial is
// not relevant and set to d(x)=x.
static POLY irred;

////////////////////////////////////////////////////////////////////////////////////////////////////
// Table of Conway polynomials for all fields (except prime fields)
// Created from GAP's Pols[] array.
////////////////////////////////////////////////////////////////////////////////////////////////////

/// @private
typedef struct {
   unsigned int q, p, n;
   POLY pol;
}
poltab_t;

static const poltab_t poltab[] = {
   {     4 , 2, 2,{1,1,1}},                            // GF(4)
   {     8 , 2, 3,{1,1,0,1}},                          // GF(8)
   {    16 , 2, 4,{1,1,0,0,1}},                        // GF(16)
   {    32 , 2, 5,{1,0,1,0,0,1}},                      // GF(32)
   {    64 , 2, 6,{1,1,0,1,1,0,1}},                    // GF(64)
   {   128 , 2, 7,{1,1,0,0,0,0,0,1}},                  // GF(128)
   {   256 , 2, 8,{1,0,1,1,1,0,0,0,1}},                // GF(256)
   {   512 , 2, 9,{1,0,0,0,1,0,0,0,0,1}},              // GF(512)
   {  1024 , 2,10,{1,1,1,1,0,1,1,0,0,0,1}},            // GF(1024)
   {  2048 , 2,11,{1,0,1,0,0,0,0,0,0,0,0,1}},          // GF(2048)
   {  4096 , 2,12,{1,1,0,1,0,1,1,1,0,0,0,0,1}},        // GF(4096)
   {  8192 , 2,13,{1,1,0,1,1,0,0,0,0,0,0,0,0,1}},      // GF(8192)
   { 17384 , 2,14,{1,0,0,1,0,1,0,1,0,0,0,0,0,0,1}},    // GF(16384)
   { 32768 , 2,15,{1,0,1,0,1,1,0,0,0,0,0,0,0,0,0,1}},  // GF(32768)
   { 65536 , 2,16,{1,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,1}},// GF(65536)
   {     9 , 3, 2,{2,2,1}},                            // GF(9)
   {    27 , 3, 3,{1,2,0,1}},                          // GF(27)
   {    81 , 3, 4,{2,0,0,2,1}},                        // GF(81)
   {   243 , 3, 5,{1,2,0,0,0,1}},                      // GF(243)
   {   729 , 3, 6,{2,2,1,0,2,0,1}},                    // GF(729)
   {  2187 , 3, 7,{1,0,2,0,0,0,0,1}},                  // GF(2187)
   {  6561 , 3, 8,{2,2,2,0,1,2,0,0,1}},                // GF(6561)
   { 19683 , 3, 9,{1,1,2,2,0,0,0,0,0,1}},              // GF(19683)
   { 59049 , 3,10,{2,1,0,0,2,2,2,0,0,0,1}},            // GF(59049)
   {    25 , 5, 2,{2,4,1}},                            // GF(25)
   {   125 , 5, 3,{3,3,0,1}},                          // GF(125)
   {   625 , 5, 4,{2,4,4,0,1}},                        // GF(625)
   {  3125 , 5, 5,{3,4,0,0,0,1}},                      // GF(3125)
   { 15625 , 5, 6,{2,0,1,4,1,0,1}},                    // GF(15625)
   {    49 , 7, 2,{3,6,1}},                            // GF(49)
   {   343 , 7, 3,{4,0,6,1}},                          // GF(343)
   {  2401 , 7, 4,{3,4,5,0,1}},                        // GF(2401)
   { 16807 , 7, 5,{4,1,0,0,0,1}},                      // GF(16807)
   {   121, 11, 2,{2,7,1}},                            // GF(121)
   {  1331, 11, 3,{9,2,0,1}},                          // GF(1331)
   { 14641, 11, 4,{2,10,8,0,1}},                       // GF(14641)
   {   169, 13, 2,{2,12,1}},                           // GF(169)
   {  2197, 13, 3,{11,2,0,1}},                         // GF(2197)
   { 28561, 13, 4,{2,12,3,0,1}},                       // GF(28561)
   {   289, 17, 2,{3,16,1}},                           // GF(289)
   {  4913, 17, 3,{14,1,0,1}},                         // GF(4913)
   {   361, 19, 2,{2,18,1}},                           // GF(361)
   {  6859, 19, 3,{17,4,0,1}},                         // GF(6859)
   {   529, 23, 2,{5,21,1}},                           // GF(529)
   { 12167, 23, 3,{18,2,0,1}},                         // GF(12167)
   {   841, 29, 2,{2,24,1}},                           // GF(841)
   { 24389, 29, 3,{27,2,0,1}},                         // GF(24389)
   {   961, 31, 2,{3,29,1}},                           // GF(961)
   { 29791, 31, 3,{28,1,0,1}},                         // GF(29791)
   {  1369, 37, 2,{2,33,1}},                           // GF(1369)
   { 50653, 37, 3,{35,6,0,1}},                         // GF(50653)
   {  1681, 41, 2,{6,38,1}},                           // GF(1681)
   {  1849, 43, 2,{3,42,1}},                           // GF(1849)
   {  2209, 47, 2,{5,45,1}},                           // GF(2209)
   {  2809, 53, 2,{2,49,1}},                           // GF(2809)
   {  3481, 59, 2,{2,58,1}},                           // GF(3481)
   {  3721, 61, 2,{2,60,1}},                           // GF(3721)
   {  4489, 67, 2,{2,63,1}},                           // GF(4489)
   {  5041, 71, 2,{7,69,1}},                           // GF(5041)
   {  5329, 73, 2,{5,70,1}},                           // GF(5329)
   {  6241, 79, 2,{3,78,1}},                           // GF(6241)
   {  6889, 83, 2,{2,82,1}},                           // GF(6889)
   {  7921, 89, 2,{3,82,1}},                           // GF(7921)
   {  9409, 97, 2,{5,96,1}},                           // GF(9409)
   {10201, 101, 2,{2,97,1}},                           // GF(10201)
   {10609, 103, 2,{5,102,1}},                          // GF(10609)
   {11449, 107, 2,{2,103,1}},                          // GF(11449)
   {11881, 109, 2,{6,108,1}},                          // GF(11881)
   {12769, 113, 2,{3,101,1}},                          // GF(12769)
   {16129, 127, 2,{3,126,1}},                          // GF(16129)
   {17161, 131, 2,{2,127,1}},                          // GF(17161)
   {18769, 137, 2,{3,131,1}},                          // GF(18769)
   {19321, 139, 2,{2,138,1}},                          // GF(19321)
   {22201, 149, 2,{2,145,1}},                          // GF(22201)
   {22801, 151, 2,{6,149,1}},                          // GF(22801)
   {24649, 157, 2,{5,152,1}},                          // GF(24649)
   {26569, 163, 2,{2,159,1}},                          // GF(26569)
   {27889, 167, 2,{5,166,1}},                          // GF(27889)
   {29929, 173, 2,{2,169,1}},                          // GF(29929)
   {32041, 179, 2,{2,172,1}},                          // GF(32041)
   {32761, 181, 2,{2,177,1}},                          // GF(32761)
   {36481, 191, 2,{19,190,1}},                         // GF(36481)
   {37249, 193, 2,{5,192,1}},                          // GF(37249)
   {38809, 197, 2,{2,192,1}},                          // GF(38809)
   {39601, 199, 2,{3,193,1}},                          // GF(39601)
   {44521, 211, 2,{2,207,1}},                          // GF(44521)
   {49729, 223, 2,{3,221,1}},                          // GF(49729)
   {51529, 227, 2,{2,220,1}},                          // GF(51529)
   {52441, 229, 2,{6,228,1}},                          // GF(52441)
   {54289, 233, 2,{3,232,1}},                          // GF(54289)
   {57121, 239, 2,{7,237,1}},                          // GF(57121)
   {58081, 241, 2,{7,238,1}},                          // GF(58081)
   {63001, 251, 2,{6,242,1}},                          // GF(63001)
   { 0, 0, 0,{0}}                                      // End of table
};

static void lookupDefiningPolynomial(poltab_t* result, int q)
{
   memset(result,0,sizeof(*result));
   int p = 2;
   if (q % 2 != 0) {
      for (p = 3; q % p != 0; p += 2);
   }
   int n = 1;
   int pn;
   for (pn = p; pn < q; ++n, pn *= p);
   MTX_ASSERT(pn == q);
   result->n = n;
   result->p = p;
   result->q = q;
   if (n == 1) {
      result->pol[0] = 0;
      result->pol[1] = 1;
   } else {
      const poltab_t* tp = poltab;
      for (tp = poltab; tp->n != 0 && (tp->p != p || tp->n != n); ++tp);
      MTX_ASSERT(tp->n != 0);
      memcpy(result, tp, sizeof(*result));
   }
}

// Stores the generating (Conway) polynomial for F(«Q») in «irred».

static void getpol()
{
   poltab_t entry;
   lookupDefiningPolynomial(&entry, Q);
   MTX_ASSERT(entry.q == Q);
   memcpy(irred, entry.pol, sizeof(irred));
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Prints a polynomial at message level 1.

static void printpol(POLY a)
{
   int i,flag = 0;

   for (i = MAXPWR; i >= 0; i--) {
      if (a[i] != 0) {
         if (flag) { MESSAGE(1,("+"));}
         if (a[i] != 1) { MESSAGE(1,("%d",(int)a[i]));}
         MESSAGE(1,("x^%d",i));
         flag = 1;
      }
   }
   MESSAGE(1,("\n"));
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Inserts p into the polynomial (calculation in Z[x]).

static int polInsert(POLY a, int p)
{
   uint32_t k = 0;
   for (int i = MAXPWR; i >= 0; i--) {
      k = k * p + a[i];
   }
   return k;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Multiplies a polynomial by X.

static void polmultx(POLY a)
{
   MTX_ASSERT(a[MAXPWR] == 0);
   for (int i = MAXPWR; i > 0; --i) {
      a[i] = a[i - 1];
   }
   a[0] = 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Reduces the polynomial a modulo b. b must be normalized.

static void polmod(POLY a, POLY b)
{
   int i, l, dl, f;

   // l= index of leading coeff. in b (must be 1)
   for (l = MAXPWR; b[l] == 0 && l > 0; l--) {
   }
   MTX_ASSERT(b[l] == 1);

   for (dl = MAXPWR; dl >= l; dl--) {
      f = (int) a[dl];
      if (f == 0) {continue;}
      f = (int)P - f;
      for (i = 0; i <= l; ++i) {
         a[i + dl - l] =
            (uint16_t) ((f * b[i] + a[i + dl - l]) % (int)P);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int intLog(int base, int n)
{
   int l = 0;
   int bl = 1;
   while (bl < n) {
      bl *= base;
      ++l;
   }
   MTX_ASSERT(bl == n);
   return l;
}

//static unsigned int intExp(unsigned int base, unsigned int n)
//{
//   unsigned int l = 1;
//   while (n > 0) {
//      l *= base;
//      --n;
//   }
//   return l;
//}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Selects the generator for prime fields (the smallest natural number which generates Zp*).

static int selectGenerator(int p)
{
   for (unsigned gen = 1; gen < P; ++gen) {
      uint32_t x = gen;
      int i;
      for (i = 1; x != 1; ++i) {
         x = (x * gen) % P;
      }
      if (i == P - 1)
         return gen;
   }
   
   mtxAbort(MTX_HERE, "Something is wrong...");
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void computeFieldMapP(uint16_t *map, int p)
{
   int gen = selectGenerator(p);
   uint16_t aExt = 1;
   for (uint16_t a = 0; a < p - 1; ++a) {
      map[a] = aExt;
      map[p + aExt] = a;
      aExt = ((uint32_t)aExt * gen) % p;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Verify that the chosen polynomial was primitive.

static void testprim()
{
   uint32_t i;
   uint32_t *a;

   a = NALLOC(uint32_t, Q);
   memset(a, 0, sizeof(uint32_t) * (size_t)Q);
   for (i = 1; i < Q; i++) {
      MTX_ASSERT(FfFromIntTable[i] < Q);
      ++a[FfFromIntTable[i]];
   }
   for (i = 0; i < Q - 1; i++) {
      if (a[i] != 1) {
         mtxAbort(MTX_HERE, "Internl error: polynomial is not primitive");
      }
   }
   free(a);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void computeFieldMapQ(uint16_t *map, int q)
{
   uint16_t* ffToExt = map;
   uint16_t* extToFf = map + q;
   poltab_t irred;

   lookupDefiningPolynomial(&irred, q);

   // map nonzero elements to 1 ... q-1
   POLY apol;
   memset(apol,0,sizeof(POLY));
   apol[0] = 1; // start with unit element

   for (FEL a = 0; a < (FEL)(q - 1); a++) {
      const int aExt = polInsert(apol, irred.p);
      MTX_ASSERT(aExt > 0 && aExt < q);
      ffToExt[a] = aExt;
      extToFf[aExt] = a;

      // multiply with generator (=x)
      polmultx(apol);
      polmod(apol,irred.pol);
   }
   testprim();
   Gen = P;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Computes the mapping between internal and external representation for F(q) with char(F) = p.
//
// The returned pointer is a table of size 2*q, where the first half defines the internal to
// external mapping, and the second half defines the external to internal mapping:
//
//  - map[a]     = external representation ("number") of the nonzero field element a.
//                 The mapping of FF_ZERO to 0 is not in the table!
//  - map[q + i] = field element corresponding to number i
//
// The table is allocated dynamically and must be freed by the caller.

static uint16_t* computeFieldMap(const int p, const int q)
{
   const int n = intLog(p, q);

   // Initialize with FF_ZERO
   uint16_t* tbl = NALLOC(uint16_t, 2 * q);
   for (size_t k = 0; k < q; ++k) tbl[k] = 0;
   for (size_t k = q; k < 2 * q; ++k) tbl[k] = FF_ZERO;
   if (n == 1) {
      computeFieldMapP(tbl, p);
   } else {
      computeFieldMapQ(tbl, q);
   }
   return tbl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

// Computes the embedding and restriction table for F(r) < F(Q).
// «table» must point to a buffer of size 2r. The embedding is written to the first half of «table»,
// the restriction is written to the second half.

static void computeEmbedding(uint16_t r, uint16_t* table)
{
   MTX_ASSERT(r % P == 0);
   const int m = intLog(P, r);  // r = P^m
   MTX_ASSERT(N % m ==  0);

   // Fill table with "invalid" marker. FF_ZERO is handled in code and never occurs as a
   // table entry.
   uint16_t* const emb = table;
   uint16_t* const restr = table + r;
   for (size_t i = 0; i < r + Q; ++i) table[i] = FF_ZERO;

   // compute embedding + restriction (for nonzero elements)
   uint16_t* const subfieldMap = computeFieldMap(P, r);
   int d = (Q - 1) / (r - 1);
   MTX_ASSERT(d * (r - 1) == Q - 1);

   int i = 0;
   for (int j = 0; j < r - 1; i += d, ++j) {
       MTX_ASSERT(j >= 0 && j < r);
       MTX_ASSERT(i >= 0 && i < Q);
       emb[j] = i;
       restr[i] = j;
   }

   free(subfieldMap);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void computeEmbeddingTables()
{
   // Determine subfield orders and calculate the embeddings table size.
   embeddingTablesSize = 0;
   numberOfSubfields = 0;
   int subfieldOrder = P;
   for (int m = 1; m < N; ++m, subfieldOrder *= P) {
      if (N % m != 0) {continue;}
      subfieldOrders[numberOfSubfields++] = subfieldOrder;
      embeddingTablesSize += subfieldOrder + Q;
   }

   // Calculate embeddings
   embeddingTables = NREALLOC(embeddingTables, FEL, embeddingTablesSize);
   size_t offset = 0;
   for (uint16_t i = 0; i < numberOfSubfields; ++i) {
      computeEmbedding(subfieldOrders[i], embeddingTables + offset);
      offset += subfieldOrders[i] + Q;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Initialize index table (for Q not prime).

static void initarith()
{
   int i, elem;

   // start with a(X) = 1
   POLY a;
   memset(a,0,sizeof(POLY));
   a[0] = 1;

   FfFromIntTable[0] = 0xFFFF;
   FfToIntTable[Q] = 0; // never used

   for (i = 0; i < (int)Q - 1; i++) {
      elem = polInsert(a, P);
      MTX_ASSERT(elem > 0 && elem < Q);
      FfFromIntTable[elem] = (uint16_t) i;
      FfToIntTable[i] = (uint16_t) elem;
      polmultx(a);
      polmod(a,irred);
   }
   testprim();
   Gen = P;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Initialize index table (for Q prime).

static void initarithP()
{
   uint16_t a, gen, x;

   // Find generator
   for (gen = 1; gen < P; ++gen) {
      x = gen;
      uint32_t i;
      for (i = 1; x != 1; ++i) {
         x = (uint16_t) (((unsigned long)x * gen) % P);
      }
      if (i == P - 1) {break;}
   }
   Gen = gen;

   // Fill the index table
   FfFromIntTable[0] = 0xFFFF;
   FfToIntTable[Q] = 0; // never used
   a = 1;
   for (unsigned i = 0; i < P - 1; i++) {
      MTX_ASSERT(a > 0 && a < Q);
      FfFromIntTable[a] = (uint16_t) i;
      FfToIntTable[i] = a;
      a = (uint16_t) (((unsigned long) a * gen) % P);
   }
   testprim();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Computes a+1 for all field elements.

static void computeIncrementTable()
{
   for (uint32_t i = 1; i <= Q - 1; ++i) {
      const uint16_t j = (uint16_t)((i % P) == P - 1 ? i + 1 - P : i + 1);
      if (j == 0) {
         Minusone = FfFromIntTable[i];
         MESSAGE(1,("MinusOne=%u(0x%04x)\n", i, FfFromIntTable[i]));
      }
      inc[FfFromIntTable[i]] = FfFromIntTable[j];
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Creates the file and writes the header.

static void writeHeader()
{
   MESSAGE(1, ("Generating arithmetic tables\n"));
   MESSAGE(1, ("ZZZ version : %u\n", (unsigned) MTX_ZZZVERSION));
   MESSAGE(1, ("Field order : %u=%u^%u\n", (unsigned)Q, (unsigned)P, (unsigned) N));
   if (P != Q) {
      MESSAGE(1, ("Polynomial  : "));
      printpol(irred);
      MESSAGE(1, ("Generator   : x\n"));
   }
   else {
      MESSAGE(1, ("Generator   : %u\n", (unsigned) Gen));
   }

   char fname[50];
   sprintf(fname, "p%5.5ld.zzz", (long) Q);
   fd = sysFopen(fname, "wb::lib");

   const uint32_t header[5] = { MTX_ZZZVERSION, P, Q, N, Gen };
   sysWrite32(fd, header, 5);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int fieldOrder)
{
   MTX_ASSERT(FF_ZERO == 0xFFFF);

   if (fieldOrder < 2 || fieldOrder > 65536) {
      mtxAbort(MTX_HERE, "Field order %d out of range (2-65536)", fieldOrder);
   }
   Q = (uint32_t) fieldOrder;
   for (P = 2; Q % P != 0; ++P);
   uint32_t q = Q;
   for (N = 0; (q % P) == 0; q /= P, ++N);
   if (q != 1) {
      mtxAbort(MTX_HERE, "Field order %d is not allowed", fieldOrder);
   }

   inc = NREALLOC(inc, uint16_t, Q);
   FfToIntTable = NREALLOC(FfToIntTable, uint16_t, Q);
   FfFromIntTable = NREALLOC(FfFromIntTable, uint16_t, Q);
   getpol();
   if (N != 1) {
      initarith();
   } else {
      initarithP();
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void writeTables()
{
   sysWrite16(fd, &Minusone, 1);
   sysWrite16(fd, inc,(size_t)Q - 1);
   sysWrite16(fd, FfToIntTable,(size_t)Q);
   sysWrite16(fd, FfFromIntTable,(size_t)Q);
   sysWrite16(fd, &numberOfSubfields,1);
   sysWrite16(fd, subfieldOrders, numberOfSubfields);
   sysWrite16(fd, embeddingTables, embeddingTablesSize);
   fclose(fd);
   MESSAGE(1,("Ok\n"));
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int ffMakeTables(int field)
{
   const int context = mtxBegin(MTX_HERE, "Generating arithmetic tables for GF(%d)", field);
   MTX_ASSERT(sizeof(int) >= 4); // required to avoid overflows
   MTX_ASSERT(sizeof(FEL) == 2);
   init(field);
   computeIncrementTable();
   computeEmbeddingTables();
   writeHeader();
   writeTables();
   mtxEnd(context);
   return 0;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
