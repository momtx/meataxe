////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Bit Strings, basic function
//
// (C) Copyright 1998-2014 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

#define BPL (sizeof(long) * 8)          /* Number of bits in a long */

/// @addtogroup bs
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Clear a bit string
/// This function clears all bits in a bit string.
/// @return 0 on success, -1 on error.

int BsClearAll(BitString_t *bs)
{
   if (!BsIsValid(bs)) {
      return -1;
   }
   memset(bs->Data,0,bs->BufSize * sizeof(long));
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Compare two bit strings
/// @param a First bit string.
/// @param b Second bit string.
/// @return 0 if the bit strings are equal. Otherwise the return value is different from zero.

int BsCompare(const BitString_t *a, const BitString_t *b)
{
   int i;
   if (!BsIsValid(a) || !BsIsValid(b)) {
      return -1;
   }
   i = a->Size - b->Size;
   if (i != 0) {
      return i;
   }
   return memcmp(a->Data,b->Data,a->BufSize * sizeof(long));
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Copy a bit string
/// @see BsDup(BitString_t const *)
/// @param dest Destination bit string.
/// @param src Source bit string.
/// @return @em dest on success, 0 on error.

BitString_t *BsCopy(BitString_t *dest, const BitString_t *src)
{
   if (!BsIsValid(dest) || !BsIsValid(src)) {
      return NULL;
   }
   if (dest->Size != src->Size) {
      MTX_ERROR1("%E",MTX_ERR_INCOMPAT);
      return NULL;
   }
   memcpy(dest->Data,src->Data,src->BufSize * sizeof(long));
   return dest;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// @fn int BsSet(BitString_t *bs, int i)
/// Set one bit in a bit string.
/// This function sets a single bit in a bit string. @em i is the index of the bit to set,
/// the first bit having index 0. @em i must be less than the bit string's size.
/// It is not possible to change the size of an existing bit string.
/// @see BsClear(BitString_t*,int)
/// @param bs Pointer to the bit string.
/// @param i Index of the bit to set (starting from 0).
/// @return 0 on success, -1 on error.

#ifdef DEBUG

int BsSet(BitString_t *bs, int i)
{
   if (!BsIsValid(bs)) {
      return -1;
   }
   if ((i < 0) || (i >= bs->Size)) {
      MTX_ERROR2("i=%d: %E",i,MTX_ERR_BADARG);
      return -1;
   }
   bs->Data[i / BPL] |= 1L << (i % BPL);
   return 0;
}


#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @fn int BsClear(BitString_t *bs, int i)
/// Clear one bit in a bit string.
/// This function clears a single bit in a bit string. @em i is the index of the bit to clear,
/// the first bit having index 0. @em i must be less than the bit string's size.
/// It is not possible to change the size of an existing bit string.
/// @see BsSet(BitString_t*,int)
/// @param bs Pointer to the bit string.
/// @param i Index of the bit to clear (starting from 0).
/// @return 0 on success, -1 on error.

#ifdef DEBUG

int BsClear(BitString_t *bs, int i)
{
   if (!BsIsValid(bs)) {
      return -1;
   }
   if ((i < 0) || (i >= bs->Size)) {
      MTX_ERROR2("i=%d: %E",i,MTX_ERR_BADARG);
      return -1;
   }
   bs->Data[i / BPL] &= ~(1L << (i % BPL));
   return 0;
}


#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @fn int BsTest(const BitString_t *, int)
/// Test a bit in a bit string.
/// @param bs Pointer to the bit string.
/// @param i Index of the bit to test (starting from 0).
/// @return 1 if the bit is set, 0 otherwise.

#ifdef DEBUG

int BsTest(const BitString_t *bs, int i)
{
   if (!BsIsValid(bs)) {
      return -1;
   }
   if ((i < 0) || (i >= bs->Size)) {
      MTX_ERROR2("i=%d: %E",i,MTX_ERR_BADARG);
      return -1;
   }
   return (bs->Data[i / BPL] & (1L << (i % BPL))) != 0 ? 1 : 0;
}


#endif

/// @}
