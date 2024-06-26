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
/// Given a finitely generated matrix algebra A, the word generator produces a sequence of
/// "random" elements of A, i.e., words in the generators.
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
/// For a given number of generators, the computation depends only on the word number.
/// For example, word 1833 in two generators a and b is always
/// b<sup>5</sup>+aba<sup>3</sup>+a<sup>2</sup>ba<sup>2</sup>+a<sup>3</sup>ba+a<sup>4</sup>.
///
/// <b>Implementation details</b><br>
/// The generator produces words in blocks of 238, i.e., words 1 to 238 belong to block 1,
/// words 239 to 476 to block 2 and so on.
/// For each block, 8 monomials A,B,...,H in the generators are chosen by calculating "random"
/// products of the generators. Then, all possible sums of 2 up to 6 of the monomials are
/// calculated, yielding 238 words.
/// The order in which these sums are taken is fixed:
/// A+B+C, A+B+C+F, A+D+E+G+H, A+B+D+E+G+H, ..., C+D+E+F+G+H.
/// See the @c BitTab[] array in wgen.c for the complete list.
///
/// Remark: since the generators are often invertible and the word generator is typically used
/// to find words with a small but nontrival kernel, it is a good idea to take at least two
/// summands. There seems to be no reason, however, why sums of 7 or 8 monomials are not used.
///
/// The calculation of A…H involves a simple pseudorandom number generator which is seeded with
/// the block number, and some magic including the use of fixed recipes for the first two blocks
/// (words 1 to 476). The number of factors in any monomial is limited to 5 for the first 200
/// blocks, to 6 for blocks 200 to 1999, and to 7 for blocks 2000 to 19999. For example, assuming
/// two generators a and b, the summand abababa has 7 factors and thus cannot appear before block
/// 2000, which means not before word 476000.

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @class WgData_t
/// @brief
/// Word Generator Data.
/// This structure is used by the word generator to store the generators and internal data.

#define MINLEN  5
#define MTX_WG_MAXLEN  8

static const int B0Tab[8][MTX_WG_MAXLEN + 1] = {
   {0,-1,-1,-1,-1}, {1,-1,-1,-1,-1}, {2,3,-1,-1,-1}, {5,4,-1,-1,-1},
   {7,9,-1,-1,-1}, {6,11,13,-1,-1}, {17,8,1,-1,-1}, {19,10,21,23,-1}
};

static const int B1Tab[8][MTX_WG_MAXLEN + 1] = {
   {0,1,2,-1,-1}, {4,3,5,6,-1},{8,10,7,-1,-1},{12,9,14,11,-1},
   {13,16,15,18,-1}, {17,20,22,19,-1},{24,21,23,25,26,-1},{27,29,28,31,33,-1}
};

static const uint8_t BitTab[238] = {
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

static int CalcLen(int blk)
{
   if (blk < 200) {
      return 5;
   }
   if (blk < 2000) {
      return 6;
   }
   if (blk < 20000) {
      return 7;
   }
   return MTX_WG_MAXLEN;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Determines the recipe for blocks 2 and higher.
/// @param blk Block number
/// @param[out] buf The recipe (entries must still be reduced modulo N,
///    the number of generators, see MakeBuf()!).

static void MakeBufX2(uint32_t blk, int buf[8][MTX_WG_MAXLEN + 1])
{
   MTX_ASSERT(blk >= 2);
   uint32_t r = blk;
   int i;
   int m[8];
   int len, mod;

   len = CalcLen(blk);
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
/// @param blk Block number
/// @param ngen Number of generators
/// @param[out] buf The recipe. buf[i] contains the list of generators
///     that must be multiplied to calculate the i-th monomial. The list
///     is terminated by -1.
static void MakeBuf2(int buf[8][MTX_WG_MAXLEN + 1], uint32_t blk, int ngen)

{
   if (blk == 0) {
      memcpy(buf,B0Tab,sizeof(B0Tab));
   } else if (blk == 1) {
      memcpy(buf,B1Tab,sizeof(B1Tab));
   } else {
      MakeBufX2(blk,buf);
   }
   
   for (int pos = 0; pos < 8; ++pos) {
      for (int* x = buf[pos]; *x >= 0; ++x) {
         *x %= ngen;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Determines the 8 monomials for a given block of words.
/// @param blk Block number
/// @param ngen Number of generators
/// @param[out] buf The recipe. buf[i] contains the list of generators
///     that must be multiplied to calculate the i-th monomial. The list
///     is terminated by -1.

static void MakeBuf(WgData_t* wg, int blk)
{
   return MakeBuf2(wg->buf, blk, wg->Rep->NGen);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Determines one of the 8 monomials for a given block of words.
/// @param blk Block number
/// @param pos Monomial to calculate (0..7).
/// @param ngen Number of generators.
/// @return Description of the monomial, a list of generator numbers terminated by -1.
///    To calcuate the monomial, multiply these generators in the given order.

static const int *B(WgData_t* wg, int blk, int pos)
{
   if (blk != wg->lastn2) {
      // Block changed. Recalculate recipe.
      wg->lastn2 = blk;
      MakeBuf(wg, blk);
   }
   return wg->buf[pos];
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void GenBasis(WgData_t *b, uint32_t blk, int pos)
{
   const int *x;
   Matrix_t *buf = NULL;

   if (b->Basis[pos] != NULL) {
      matFree(b->Basis[pos]);
   }
   for (x = B(b, blk, pos); *x >= 0; ++x) {
      MTX_ASSERT(*x >= 0 && *x < b->Rep->NGen);
      if (buf == NULL) {
         buf = matDup(b->Rep->Gen[*x]);
      } else {
         matMul(buf,b->Rep->Gen[*x]);
      }
   }
   MTX_ASSERT(buf != NULL);
   b->Basis[pos] = buf;
   b->N2[pos] = blk;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static Matrix_t* makeMonomial(const WgData_t *wg, uint32_t blk, int pos)
{
   Matrix_t *m = NULL;
   MTX_ASSERT(pos >= 0 && pos < 8);

   int buf[8][MTX_WG_MAXLEN + 1];
   MakeBuf2(buf, blk, wg->Rep->NGen);
   for (const int* x = buf[pos]; *x >= 0; ++x) {
      MTX_ASSERT(*x >= 0 && *x < wg->Rep->NGen);
      if (m == NULL) {
         m = matDup(wg->Rep->Gen[*x]);
      } else {
         matMul(m,wg->Rep->Gen[*x]);
      }
   }
   return m;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Symbolic name of a word.
/// This function returns a symbolic representation of the word @p n as a polynomial in the
/// generators. For example, word 306 with two generators is represented as "ab2a+a2b+ab3a".
/// The return value is a pointer to an internal buffer in the word generator, which is
/// overwritten on each call for the same word generator.
///
/// See also @ref wgDescribeWord.
///
/// @note This function is not threadsafe!
///
/// @param wg Pointer to word generator data.
/// @param n Word number.
/// @return Symbolic name of the word.

const char *wgSymbolicName(WgData_t *wg, long n)
{
   MTX_ASSERT(n > 0);

   char *c = wg->name;
   wgDescribeWord(wg,n);
   for (int* x = wg->Description; *x != -1; ) {
      if (x != wg->Description) { *c++ = '+'; }
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
   return wg->name;
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
/// @param blk Block number
/// @param i Monomial number (0..7)

static void DescribeMonomial(WgData_t *wg, int *pos, long blk, int i)
{
   const int *x = B(wg, blk, i);
   while (*x >= 0) {
      AppendDescription(wg,pos,*x++);    // Generator number
   }
   AppendDescription(wg,pos,-1);         // End of monomial
}

static void splitWordNumber(uint32_t n, uint8_t *mask, uint32_t* blk)
{
   MTX_ASSERT(n > 0);
   --n;
   *mask = BitTab[n % 238];
   *blk = n / 238;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Creates a symbolic description of a word.
/// Stores the description of the given word in «wg->Description». The
/// description is a sequence of monomials terminated by -1. Each monomial
/// itself is a sequence of integers, again terminated by -1, specifying the
/// generators that must be multiplied to obtain the monomial. The word is
/// the sum of the monomials.
///
/// For example a+b+baa would be represented as 0,-1,1,-1,1,0,0,-1,-1.
///
/// See also @ref wgSymbolicName.
///
/// @note This function is not threadsafe!
/// «wg->Description» is overwritten each time this function is called and
/// should treated as read-only. In particular, do not attempt to free or
/// reallocate the buffer!
///
/// @param wg Word generator.
/// @param n Word number (> 0).

int *wgDescribeWord(WgData_t *wg, uint32_t n)
{
   int pos = 0;

   uint8_t n1;
   uint32_t blk;
   splitWordNumber(n, &n1, &blk);

   for (int i = 0; i < 8 && n1 != 0; ++i, n1 >>= 1) {
      if ((n1 % 2) != 0) {
         DescribeMonomial(wg,&pos,blk,i);
      }
   }
   AppendDescription(wg,&pos,-1); // end of description
   return wg->Description;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Calculates a word.
/// This function calculates an element in the algebra generated by a set of matrices. The element
/// is identified by a single number. If WgData_t structures were initialized with representations
/// of the same group generators, both generators will produce representations of the same group
/// algebra element for any fixed number.
///
/// @note This function is not threadsafe. Use @ref wgMakeWord2 to use the same word generator
/// from multiple threads.
///
/// @param wg Pointer to word generator data.
/// @param n Word number.
/// @return Matrix representation of the word

Matrix_t *wgMakeWord(WgData_t *wg, uint32_t n)
{
   MTX_ASSERT(n > 0);
   uint8_t n1;
   uint32_t blk;
   splitWordNumber(n, &n1, &blk);

   Matrix_t *w = NULL;
   for (int i = 0; i < 8 && n1 != 0; ++i, n1 >>= 1) {
      if (n1 % 2 == 0) {
         continue;
      }
      if (wg->N2[i] != blk) {
         GenBasis(wg,blk,i);
      }
      if (w == NULL) {
         w = matDup(wg->Basis[i]);
      } else {
         matAdd(w,wg->Basis[i]);
      }
   }
   return w;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Calculates a word (threadsafe version)
/// This function works like wgMakeWord() but does not access the internal state of the word
/// generator. It may be used in different threads with the same WgData_t structure.

Matrix_t *wgMakeWord2(WgData_t *wg, uint32_t n)
{
   MTX_ASSERT(n > 0);
   uint8_t mask;
   uint32_t blk;
   splitWordNumber(n, &mask, &blk);
   
   Matrix_t *w = NULL;
   for (int i = 0; i < 8 && mask != 0; ++i, mask >>= 1) {
      if (mask % 2 == 0) {
         continue;
      }
      Matrix_t* m = makeMonomial(wg, blk, i);
      if (w == NULL) {
         w = m;
      } else {
         matAdd(w,m);
         matFree(m);
      }
   }
   return w;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void wgValidate(const struct MtxSourceLocation* where, WgData_t* wg)
{
   if (wg == NULL) {
      mtxAbort(where, "Invalid word generator (NULL)");
   }
   mrValidate(where, wg->Rep);
   for (int i = 0; i < 8; ++i) {
      if (wg->Basis[i] != 0) {
         matValidate(where, wg->Basis[i]);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Creates a word generator for a given matrix representation.
///
/// There must be at least one generator in @p rep.
///
/// The work generator only stores a reference to the generators but not take ownership.
/// The generators must not be modified or destroyed while the word generator is alive.
/// The caller remains responsible for destroying the generators after(!) the word generator
/// was destroyed.

WgData_t *wgAlloc(const MatRep_t *rep)
{
   mrValidate(MTX_HERE, rep);

   WgData_t *wg = (WgData_t*) mmAlloc(MTX_TYPE_WORD_GENERATOR, sizeof(WgData_t));
   wg->Rep = rep;
   for (int k = 0; k < 8; ++k) {
      wg->Basis[k] = NULL;
      wg->N2[k] = -1;
   }
   wg->lastn2 = -1;
   return wg;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Destroys a word generator and releases internal resources.
/// Note: the matrix representation of the generators is not released.
/// See also @ref wgAlloc.

int wgFree(WgData_t *wg)
{
   wgValidate(MTX_HERE, wg);

   for (int k = 0; k < 8; ++k) {
      if (wg->Basis[k] != NULL) {
         matFree(wg->Basis[k]);
      }
   }
   if (wg->Description != 0) {
      sysFree(wg->Description - 1);
   }
   wg->Rep = NULL;
   mmFree(wg, MTX_TYPE_WORD_GENERATOR);
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Calculate finger print.
/// This function calculates the "finger print" of a module, i.e. the nullities
/// of the first 6 words.
///
/// @note This function is not threadsafe!
///
/// @param wg Word generator data.
/// @param fp Buffer for the finger print (6 numbers).

void wgMakeFingerPrint(WgData_t *wg, uint32_t fp[6])
{
   int i;
   for (i = 1; i <= 6; ++i) {
      fp[i - 1] = matNullity__(wgMakeWord(wg,i));
   }
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
