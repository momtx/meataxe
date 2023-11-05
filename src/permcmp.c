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
/// This function compares two permutations. If the permutations are equal, the return value is 0.
/// Otherwise the return value is positive if @em a is "greater" than  @em b, and negative if @em a
/// is "less" than @em b. The ordering of permutations is defined as follows:
/// - If the permutations have different degrees, the permutations with the smaller degree is
///   smaller.
/// - Otherwise, the order is defined by the lexicograhical order of the sequences
///   (a(0), … , a(n-1)) and (b(0),…,b(n-1)), where n is the degree.

int permCompare(const Perm_t *a, const Perm_t *b)
{
   // check arguments
   permValidate(MTX_HERE,a);
   permValidate(MTX_HERE,b);

   // compare degrees
   if (a->degree > b->degree) { return 1; }
   if (a->degree < b->degree) { return -1; }
  
   // compare the entries
   uint32_t* pa = a->data;
   uint32_t* pb = b->data;
   while (pa < a->data + a->degree) {
      if (*pa > *pb) return 1;
      if (*pa < *pb) return -1;
      ++pa;
      ++pb;
   }
   return 0;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
