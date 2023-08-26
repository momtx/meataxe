////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Bit Strings, count matching bits
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


static int BitCount[256] =  /* Number of '1' bits in binary representation */
{
   0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
   1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
   1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
   1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
   3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
   1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
   3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
   3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
   3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
   4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
};

/// @addtogroup bs
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Intersection count.
/// This function calculates the cardinality of the intersection of two bit strings, i.e.,
/// the number of bits that are set in both @em a and @em b. The arguments must be bit
/// strings of the same size.
/// @return Number of bits in the intersection of @em a and @em b, or -1 on error.

int bsIntersectionCount(const BitString_t *a, const BitString_t *b)
{
   // check arguments
   bsValidate(MTX_HERE, a);
   bsValidate(MTX_HERE, b);

   register const unsigned long *ap = (const unsigned long *) a->Data;
   register const unsigned long *bp = (const unsigned long *) b->Data;
   int count = 0;
   register int i;
   for (i = a->BufSize; i > 0; --i) {
      register unsigned long a = *ap++ & *bp++;
      while (a != 0) {
         count += BitCount[a & 0xFF];
         a >>= 8;
      }
   }
   return count;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
