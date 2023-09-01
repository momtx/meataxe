////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Word generator.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
//   Local data


/// @defgroup wgen The Word Generator
/// @{
///
/// @details
/// Given a finitely generated matrix algebra A, the word generator produces
/// a sequence of "random" elements of A, i.e., words in the generators.
/// Words are numbered starting with 1.
/// Here is an example demonstrating the usage of the word generator:
/// @code
/// MatRep_t *rep;
/// ...
/// WgData_t *wg = wgAlloc(rep);
/// Matrix_t *word = wgMakeWord(wg,1833);
/// long nul = matNullity__(word);
/// printf("Word 1833 has nullity %ld\n",nul);
/// wgFree(wg);
/// @endcode
///
/// The generator produces words in blocks of 238, i.e.,
/// words 1 to 238 belong to block 1, words 239 to 476 to block 2 and so on.
/// For each block, 8 monomials a,b,...,h in the generators are chosen by
/// calculating "random" products of the generators. Then, all possible
/// sums of 2 up to 6 of the monomials are calculated, yielding 238 words.
/// The order in which these sums are taken is fixed:
/// a+b+c, a+b+c+f, a+d+e+g+h, a+b+d+e+g+h, ..., c+d+e+f+g+h.
/// See the @c BitTab[] array in wgen.c for the complete list.
///
/// Remark: since the generators are often invertible and the word generator is
/// typically used to find words with a small but nontrival kernel,
/// it is a good idea to take at least two summands. There seems to be no
/// reason, however, why sums of 7 or 8 monomials are not used.
///
/// The calculation of a...h involves a simple pseudorandom number generator
/// which is seeded with the block number, and some magic including the use
/// of fixed recipes for the first two blocks (words 1 to 476).
/// The number of factors in any monomial is limited to 5 for
/// the first 200 blocks, to 6 for blocks 200 to 1999, and to 7 for
/// blocks 2000 to 19999. For example, assuming two generators x and y,
/// the summand xyxyxyx has 7 factors and thus cannot appear before
/// block 2000, which means not before word 476000.

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @class WgData_t
/// @brief
/// Word Generator Data.
/// This structure is used by the word generator to store the generators and
/// internal data.

#define MINLEN  5
#define MAXLEN  8

static int B0Tab[8][MAXLEN + 1] = {
   {0,-1,-1,-1,-1}, {1,-1,-1,-1,-1}, {2,3,-1,-1,-1}, {5,4,-1,-1,-1},
   {7,9,-1,-1,-1}, {6,11,13,-1,-1}, {17,8,1,-1,-1}, {19,10,21,23,-1}
};

static int B1Tab[8][MAXLEN + 1] = {
   {0,1,2,-1,-1}, {4,3,5,6,-1},{8,10,7,-1,-1},{12,9,14,11,-1},
   {13,16,15,18,-1}, {17,20,22,19,-1},{24,21,23,25,26,-1},{27,29,28,31,33,-1}
};

static int BitTab[238] = {
   0x07,0x27,0xD9,0xDB,0xDF,0xF9,0xE0,0x03,0x05,0x06,0x09,0x0A,0x0B,0x0C,
   0x0D,0x0E,0x0F,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x18,0x19,0x1A,0x1B,
   0x1C,0x1D,0x1E,0x1F,0x21,0x22,0x23,0x24,0x25,0x26,0x28,0x29,0x2A,0x2B,
   0x2C,0x2D,0x2E,0x2F,0x30,0x31,0x32,0x33,0x34,0x35,0x36,0x37,0x38,0x39,
   0x3A,0x3B,0x3C,0x3D,0x3E,0x3F,0x41,0x42,0x43,0x44,0x45,0x46,0x47,0x48,
   0x49,0x4A,0x4B,0x4C,0x4D,0x4E,0x4F,0x50,0x51,0x52,0x53,0x54,0x55,0x56,
   0x57,0x58,0x59,0x5A,0x5B,0x5C,0x5D,0x5E,0x5F,0x60,0x61,0x62,0x63,0x64,
   0x65,0x66,0x67,0x68,0x69,0x6A,0x6B,0x6C,0x6D,0x6E,0x6F,0x70,0x71,0x72,
   0x73,0x74,0x75,0x76,0x77,0x78,0x79,0x7A,0x7B,0x7C,0x7D,0x7E,0x81,0x82,
   0x83,0x84,0x85,0x86,0x87,0x88,0x89,0x8A,0x8B,0x8C,0x8D,0x8E,0x8F,0x90,
   0x91,0x92,0x93,0x94,0x95,0x96,0x97,0x98,0x99,0x9A,0x9B,0x9C,0x9D,0x9E,
   0x9F,0xA0,0xA1,0xA2,0xA3,0xA4,0xA5,0xA6,0xA7,0xA8,0xA9,0xAA,0xAB,0xAC,
   0xAD,0xAE,0xAF,0xB0,0xB1,0xB2,0xB3,0xB4,0xB5,0xB6,0xB7,0xB8,0xB9,0xBA,
   0xBB,0xBC,0xBD,0xBE,0xC0,0xC1,0xC2,0xC3,0xC4,0xC5,0xC6,0xC7,0xC8,0xC9,
   0xCA,0xCB,0xCC,0xCD,0xCE,0xCF,0xD0,0xD1,0xD2,0xD3,0xD4,0xD5,0xD6,0xD7,
   0xD8,0xDA,0xDC,0xDD,0xDE,0xE1,0xE2,0xE3,0xE4,0xE5,0xE6,0xE7,0xE8,0xE9,
   0xEA,0xEB,0xEC,0xED,0xEE,0xF0,0xF1,0xF2,0xF3,0xF4,0xF6,0xF8,0xFA,0xFC
};

#define RND(x) (((x) * 214013L + 2531011L) & 0xFFFFFFFF)

////////////////////////////////////////////////////////////////////////////////////////////////////

static int CalcLen(int n2)
{
   if (n2 < 200) {
      return 5;
   }
   if (n2 < 2000) {
      return 6;
   }
   if (n2 < 20000) {
      return 7;
   }
   return MAXLEN;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Determines the recipe for blocks 2 and higher.
/// @param n2 Block number
/// @param[out] buf The recipe (entries must still be reduced modulo N,
///    the number of generators, see MakeBuf()!).

static void MakeBufX2(int n2, int buf[8][MAXLEN + 1])
{
   unsigned long r = n2;
   int i;
   int m[8];
   int len, mod;

   len = CalcLen(n2);
   mod = (1 << (len + 1)) - 2;

   for (i = 0; i < 8; ) {
      int mm, k;
      r = RND(r);
      mm = (r >> 16) % mod + 2;
      for (k = 0; k < i && m[k] != mm; ++k) {
      }
      if (k >= i) {
         m[i++] = mm;
      }
   }
   for (i = 0; i < 8; ++i) {
      unsigned long x = m[i];
      int k = len;
      int mask = 1 << len;
      while ((x & mask) == 0) {
         buf[i][k] = -1;
         mask >>= 1;
         --k;
      }
      buf[i][k] = -1;
      mask >>= 1;
      --k;
      while (k >= 0) {
         int gg = (int)((r >> 16) * 2);
         r = RND(r);
         buf[i][k] = ((x & mask) == 0) ? gg : gg + 1;
         mask >>= 1;
         --k;
      }
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Determines the 8 monomials for a given block of words.
/// @param n2 Block number
/// @param ngen Number of generators
/// @param[out] buf The recipe. buf[i] contains the list of generators
///     that must be multiplied to calculate the i-th monomial. The list
///     is terminated by -1.

static void MakeBuf(int n2, int ngen, int buf[8][MAXLEN + 1])
{
   int pos;

   if (n2 == 0) {
      memcpy(buf,B0Tab,sizeof(B0Tab));
   } else if (n2 == 1) {
      memcpy(buf,B1Tab,sizeof(B0Tab));
   } else {
      MakeBufX2(n2,buf);
   }
   for (pos = 0; pos < 8; ++pos) {
      int *x;
      for (x = buf[pos]; *x >= 0; ++x) {
         *x %= ngen;
      }
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Determines one of the 8 monomials for a given block of words.
/// @param n2 Block number
/// @param pos Monomial to calculate (0..7).
/// @param ngen Number of generators.
/// @return Description of the monomial, a list of generator numbers terminated by -1.
///    To calcuate the monomial, multiply these generators in the given order.

static const int *B(int n2, int pos, int ngen)
{
   static int buf[8][MAXLEN + 1];  /* Local cache to avoid multiple calls of MakeBuf() */
   static int lastn2 = -1;

   if (n2 != lastn2) {
      /* We have a new block. Recalculate and store the recipe. */
      MakeBuf(n2,ngen,buf);
      lastn2 = n2;
   }
   return buf[pos];
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void GenBasis(WgData_t *b, int n2, int pos)
{
   const int *x;
   Matrix_t *buf = NULL;

   if (b->Basis[pos] != NULL) {
      matFree(b->Basis[pos]);
   }
   for (x = B(n2,pos,b->Rep->NGen); *x >= 0; ++x) {
      MTX_ASSERT(*x >= 0 && *x < b->Rep->NGen);
      if (buf == NULL) {
         buf = matDup(b->Rep->Gen[*x]);
      } else {
         matMul(buf,b->Rep->Gen[*x]);
      }
   }
   MTX_ASSERT(buf != NULL);
   b->Basis[pos] = buf;
   b->N2[pos] = n2;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Symbolic name of a word.
/// This function returns a symbolic representation of the word @a n as
/// a polynomial in the generators. The generators are named a, b, c...
/// The return value is a pointer to a static buffer which is overwritten on
/// each call.
/// @param b Pointer to word generator data.
/// @param n Word number.
/// @return Symbolic name of the word.

const char *wgSymbolicName(WgData_t *b, long n)
{
   static char name[8 * (MAXLEN + 1) + 1];
   char *c = name;

   MTX_ASSERT(n > 0);

   wgDescribeWord(b,n);
   int *x;
   for (x = b->Description; *x != -1; ) {
      if (x != b->Description) { *c++ = '+'; }
      do {
         int *gen = x;
         *c++ = *gen + 'a';
         while (*x == *gen) {
            ++x;
         }
         if (x - gen > 1) { *c++ = x - gen + '0'; }
      } while (*x != -1);
      ++x;
   }
   *c = 0;
   return name;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void AppendDescription(WgData_t *b, int *pos, int x)
{
   int capacity = b->Description ? b->Description[-1] : 0;
   if (*pos >= capacity) {
      capacity += 32;
      const size_t size = capacity * sizeof(int) + 1;
      b->Description = b->Description == 0 ?
                       (int*) sysMalloc(size) + 1 :
                       (int*) sysRealloc(b->Description - 1, size) + 1;
      b->Description[-1] = capacity;
   }
   b->Description[*pos] = x;
   ++*pos;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// @param b Word generator.
/// @param pos Length of the description.
/// @param n2 Block number
/// @param i Monomial number (0..7)

static void DescribeMonomial(WgData_t *b, int *pos, long n2, int i)
{
   const int ngen = b->Rep->NGen;
   const int *x = B(n2,i,ngen);
   while (*x >= 0) {
      AppendDescription(b,pos,*x++);    /* Generator number */
   }
   AppendDescription(b,pos,-1);         /* End of monomial */
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Creates a symbolic description of a word.
/// Stores the description of the given word in «b->Description». The
/// description is a sequence of monomials terminated by -1. Each monomial
/// itself is a sequence of integers, again terminated by -1, specifying the
/// generators that must be multiplied to obtain the monomial. The word is
/// the sum of the monomials.
///
/// For example a+b+baa would be represented as 0,-1,1,-1,1,0,0,-1,-1.
///
/// «b->Description» is overwritten each time this function is called and
/// should treated as read-only. In particular, do not attempt to free or
/// reallocate the buffer!
///
/// @param b Word generator.
/// @param n Wird number.

int *wgDescribeWord(WgData_t *b, long n)
{
   int i;
   int pos = 0;
   long n1, n2;

   MTX_ASSERT(n > 0);
   --n;
   n1 = BitTab[(int)(n % 238)];
   n2 = (int)(n / 238);

   for (i = 0; i < 8 && n1 != 0; ++i, n1 >>= 1) {
      if ((n1 % 2) != 0) {
         DescribeMonomial(b,&pos,n2,i);
      }
   }
   AppendDescription(b,&pos,-1);        /* End marker */
   return b->Description;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Calculates a word.
/// This function calculates an element in a matrix algebra, given its number.
/// Generators for the algebra are specified when the WgData_t structure is
/// allocated with wgAlloc().
/// @param b Pointer to word generator data.
/// @param n Word number.
/// @return Matrix representation of the word or 0 on error.

Matrix_t *wgMakeWord(WgData_t *b, long n)
{
   Matrix_t *w = NULL;
   int n1, n2;
   int i;

   MTX_ASSERT(n > 0);
   --n;
   n1 = BitTab[(int)(n % 238)];
   n2 = (int)(n / 238);
   for (i = 0; i < 8 && n1 != 0; ++i, n1 >>= 1) {
      if ((n1 % 2) == 0) {
         continue;
      }
      if (b->N2[i] != n2) {
         GenBasis(b,n2,i);
      }
      if (w == NULL) {
         w = matDup(b->Basis[i]);
      } else {
         matAdd(w,b->Basis[i]);
      }
   }
   return w;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static int CheckArgs(const MatRep_t *rep)
{
   if (!mrIsValid(rep)) {
      mtxAbort(MTX_HERE,"rep: %s",MTX_ERR_BADARG);
      return -1;
   }
   if (rep->NGen < 1) {
      mtxAbort(MTX_HERE,"Invalid number of generators (%d)",rep->NGen);
      return -1;
   }
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Initialize the word generator.
/// This function initializes the word generator for a given matrix
/// representation @a rep. There must be at least one generator.
/// On success, wgAlloc() returns a pointer to an internal data
/// structure which can be used in subsequent calls to wgMakeWord()
/// and wgFree(). If an error occurs, the return value is 0.
///
/// Note that the word generator does not create internal copies of the
/// generators. The caller must assure that the generators are not deleted
/// or modified as long as the word generator is in use.
///
/// @param rep Pointer to the matrix representation.
/// @return Pointer to a new word generator data structure or 0 on error.

WgData_t *wgAlloc(const MatRep_t *rep)
{
   int k;
   WgData_t *b;

   /* Check the Arguments
      ------------------- */
   if (CheckArgs(rep) != 0) {
      return NULL;
   }

   /* Create a new data structure
      ---------------------------  */
   if ((b = ALLOC(WgData_t)) == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate word generator data: %S");
      return NULL;
   }

   /* Initialize the members
      ---------------------- */
   b->Rep = rep;
   for (k = 0; k < 8; ++k) {
      b->Basis[k] = NULL;
      b->N2[k] = -1;
   }
   b->Description = NULL;

   return b;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Terminate the word generator.
/// @param b Pointer to word generator data.
/// @return 0 on success, -1 on error.
/// This function terminates the word generator and cleans up internal data structures.
/// Note that the generators that were passed to wgAlloc() are not freed.

int wgFree(WgData_t *b)
{
   int k;

   /* Check the handle
      ---------------- */
   if (b == NULL) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_BADARG);
      return -1;
   }

   /* Free all basis matrices
      ----------------------- */
   for (k = 0; k < 8; ++k) {
      if (b->Basis[k] != NULL) {
         matFree(b->Basis[k]);
      }
   }
   if (b->Description != 0) {
      sysFree(b->Description - 1);
   }
   memset(b,0,sizeof(WgData_t));
   sysFree(b);
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Calculate finger print.
/// This function calculates the "finger print" of a module, i.e. the nullities
/// of the first 6 words.
/// @param b Word generator data.
/// @param fp Buffer for the finger print (6 numbers).

void wgMakeFingerPrint(WgData_t *b, int fp[6])
{
   int i;
   for (i = 1; i <= 6; ++i) {
      fp[i - 1] = matNullity__(wgMakeWord(b,i));
   }
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
