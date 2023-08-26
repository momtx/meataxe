////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Bit Strings, core functions
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

/// @defgroup bs Bit Strings
/// @{
/// @details
/// The BitString_t type represents a string of 0's and 1's. You may also think of it as a subset
/// {0,1,2,...,N-1}, where each 1 in the bit string means that the subset contains the
/// corresponding element. Unlike the set_t type bit strings are static objects which cannot
/// change their size. However, operations like the intersection of two sets are much faster with
/// BitString_t than with set_t objects.
///
/// Bit strings can be written to and read from files. The file format consists
/// of a header followed by the data. The header contains three 32-bit integers.
/// The first number is $-3$, the second number is the size of the bit
/// string, and the third number is always zero. The file header is followed by
/// the bit string data which is written as a sequence of 32-bit integers where
/// bit 0 of the first integer contains the first bit of the string.

/// @class BitString_t
/// A bit string.

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


#define BS_MAGIC 0x3ff92541

#define BS(size) (((size) + sizeof(long) - 1) / sizeof(long))

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns 1 if the given bit string is valid and 0 otherwise.

int bsIsValid(const BitString_t *bs)
{
   return bs != NULL && bs->Magic == BS_MAGIC && bs->Size >= 0
         && bs->BufSize == (int) BS(bs->Size);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Checks if a bit string is valid and aborts the progam if the test fails.

void bsValidate(const struct MtxSourceLocation* src, const BitString_t *bs)
{
   if (bs == NULL) {
      mtxAbort(MTX_HERE,"NULL bit string");
   }
   if ((bs->Magic != BS_MAGIC) || (bs->Size < 0)) {
      mtxAbort(MTX_HERE,"Invalid bit string (magic=%d, size=%d)",
                 (int)bs->Magic,bs->Size);
   }
   if (bs->BufSize != (int) BS(bs->Size)) {
      mtxAbort(MTX_HERE,"Inconsistent bit string size %d, %d)",
                 bs->Size,(int) BS(bs->Size));
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Create a bit string.
/// This function creates a new bit string of the specified size. The new bit string is
/// filled with zeros.
/// @see bsFree()
/// @param size Size of the bit string.
/// @return Pointer to the new bit string, or  NULL on error.

BitString_t *bsAlloc(int size)
{
   BitString_t *n;
   int bufsize;         /* Number of long integers in data buffer */

   if (size < 0) {
      mtxAbort(MTX_HERE,"Illegal size %d",size);
      return NULL;
   }

   // allocate data buffer
   bufsize = BS(size);
   n = (BitString_t *) sysMalloc(sizeof(BitString_t) +
                                 (bufsize == 0 ? 0 : bufsize - 1) * sizeof(long));
   if (n == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate bit string");
      return NULL;
   }

   // initialize
   n->Magic = BS_MAGIC;
   n->Size = size;
   n->BufSize = bufsize;
   memset(n->Data,0,sizeof(long) * bufsize);

   return n;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Free a bit string.
/// This function frees a bit string, releasing all associated memory.
/// @see bsAlloc()

int bsFree(BitString_t *bs)
{
   bsValidate(MTX_HERE, bs);
   memset(bs,0,sizeof(BitString_t));
   sysFree(bs);
   return 0;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
