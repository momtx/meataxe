////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Permutation order
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup perm
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Order of a permutation.
/// @param perm Pointer to the permutation.
/// @return The order of 0, or -1 on error.

int permOrder(const Perm_t *perm)
{
   long ord = 1, tord;
   int deg;
   long *p;
   long *end;
   register long *seed, *x;

   // Check arguments
   permValidate(MTX_HERE, perm);

   if (perm->Degree < 2) {
      return 1;
   }
   deg = perm->Degree;
   p = perm->Data;
   end = p + deg;
   seed = p - 1;

   // Calculate the order by running through all orbits. To keep track of
   // which points have been visited, we flip the sign of the corresponding
   // entries in the permuatation. With this trick we don't need any
   // variable-sized buffer, but we must undo the changes afterwards.
   while (1) {
      // Find next starting point
      while (*++seed < 0 && seed != end) {
      }
      if (seed == end) {
         break;             // Done!
      }
      // Find orbit size
      tord = 0;
      x = seed;
      while (1) {
         register long tmp = *x;
         if (tmp < 0) {
            break;          // Done!
         }
         *x = -tmp - 1;
         x = p + tmp;
         ++tord;
      }
      MTX_ASSERT(tord > 0, 0);

      // Calculate the l.c.m of all orbit sizes
      ord = ord / gcd(ord,tord) * tord;
   }

   // Now, all entries have their sign changed. Undo this.
   for (x = p; x != end; ++x) {
      *x = -*x - 1;
   }

   return ord;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
