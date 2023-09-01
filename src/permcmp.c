////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Compare permutations
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup perm
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Compare two permutations.
/// This function compares two permutations. If the permutations are equal,
/// the return value is 0. Otherwise the return value is positive if @em a
/// is "greater" than  @em b, and negative if @em a is "less" than @em b. The
/// ordering of permutations is defined as follows. If the permutations have
/// different degrees, the permutations with the smaller degree is smaller.
/// Otherwise, the result of the comparison is unspecified.
///
/// Note that the return value -1 does not necessarily mean that an error occured.
/// @param a Pointer to the first permutation.
/// @param b Pointer to the second permutation.
/// @return 0, if the matrices are equal, a nonzero value otherwise, -1 on error.

int permCompare(const Perm_t *a, const Perm_t *b)
{
   int i;

   // check arguments
   permValidate(MTX_HERE,a);
   permValidate(MTX_HERE,b);

   // compare degrees
   // TODO: do not return -1
   if ((i = a->Degree - b->Degree) != 0) {
      return i;
   }

   // compare the entries
   // TODO: do not return -1
   i = memcmp(a->Data,b->Data,sizeof(a->Data[0]) * a->Degree);
   if (i < 0) {
      return -1;
   } else if (i > 0) {
      return 1;
   }
   return 0;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
