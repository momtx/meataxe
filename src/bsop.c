////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Bit Strings, basic function
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


#define BPL (sizeof(long) * 8)          /* Number of bits in a long */

/// @addtogroup bs
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Clear a bit string
/// This function clears all bits in a bit string.
/// @return 0 on success, -1 on error.

int bsClearAll(BitString_t *bs)
{
   bsValidate(MTX_HERE, bs);
   memset(bs->Data,0,bs->BufSize * sizeof(long));
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Compare two bit strings
/// @param a First bit string.
/// @param b Second bit string.
/// @return 0 if the bit strings are equal. Otherwise the return value is different from zero.

int bsCompare(const BitString_t *a, const BitString_t *b)
{
   int i;
   bsValidate(MTX_HERE, a);
   bsValidate(MTX_HERE, b);
   i = a->Size - b->Size;
   if (i != 0) {
      return i;
   }
   return memcmp(a->Data,b->Data,a->BufSize * sizeof(long));
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Copy a bit string
/// @see bsDup(BitString_t const *)
/// @param dest Destination bit string.
/// @param src Source bit string.
/// @return @em dest on success, 0 on error.

BitString_t *bsCopy(BitString_t *dest, const BitString_t *src)
{
   bsValidate(MTX_HERE, src);
   bsValidate(MTX_HERE, dest);
   if (dest->Size != src->Size) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
      return NULL;
   }
   memcpy(dest->Data,src->Data,src->BufSize * sizeof(long));
   return dest;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// @fn int bsSet(BitString_t *bs, int i)
/// Set one bit in a bit string.
/// This function sets a single bit in a bit string. @em i is the index of the bit to set,
/// the first bit having index 0. @em i must be less than the bit string's size.
/// It is not possible to change the size of an existing bit string.
/// @see bsClear(BitString_t*,int)
/// @param bs Pointer to the bit string.
/// @param i Index of the bit to set (starting from 0).
/// @return 0 on success, -1 on error.

#ifdef MTX_DEBUG

int bsSet(BitString_t* bs, int i)
{
   bsValidate(MTX_HERE, bs);
   if ((i < 0) || (i >= bs->Size)) {
      mtxAbort(MTX_HERE,"i=%d: %s",i,MTX_ERR_BADARG);
      return -1;
   }
   bs->Data[i / BPL] |= 1L << (i % BPL);
   return 0;
}

#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @fn int bsClear(BitString_t *bs, int i)
/// Clear one bit in a bit string.
/// This function clears a single bit in a bit string. @em i is the index of the bit to clear,
/// the first bit having index 0. @em i must be less than the bit string's size.
/// It is not possible to change the size of an existing bit string.
/// @see bsSet(BitString_t*,int)
/// @param bs Pointer to the bit string.
/// @param i Index of the bit to clear (starting from 0).
/// @return 0 on success, -1 on error.

#ifdef MTX_DEBUG

int bsClear(BitString_t *bs, int i)
{
   bsValidate(MTX_HERE, bs);
   if ((i < 0) || (i >= bs->Size)) {
      mtxAbort(MTX_HERE,"i=%d: %s",i,MTX_ERR_BADARG);
      return -1;
   }
   bs->Data[i / BPL] &= ~(1L << (i % BPL));
   return 0;
}


#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @fn int bsTest(const BitString_t *, int)
/// Test a bit in a bit string.
/// @param bs Pointer to the bit string.
/// @param i Index of the bit to test (starting from 0).
/// @return 1 if the bit is set, 0 otherwise.

#ifdef MTX_DEBUG

int bsTest(const BitString_t *bs, int i)
{
   bsValidate(MTX_HERE, bs);
   if ((i < 0) || (i >= bs->Size)) {
      mtxAbort(MTX_HERE,"i=%d: %s",i,MTX_ERR_BADARG);
      return -1;
   }
   return (bs->Data[i / BPL] & (1L << (i % BPL))) != 0 ? 1 : 0;
}


#endif

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
