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

uint32_t permOrder(const Perm_t *perm)
{
   uint32_t order = 1;

   permValidate(MTX_HERE, perm);
   if (perm->Degree < 2) {
      return 1;
   }

   const uint32_t deg = perm->Degree;
   const uint32_t* p = perm->Data;

   uint8_t* done = NALLOC(uint8_t, deg);
   uint8_t* seed = done - 1;
   uint8_t* seedEnd = done + deg;

   // Calculate the order by running through all orbits. 
   while (1) {
      // Find start of next orbit.
      while (*++seed != 0 && seed != seedEnd);
      if (seed == seedEnd)
         break;             // Done!
      
      // Find orbit size
      uint32_t orbitSize = 0;
      uint32_t x = (uint32_t)(seed - done);
      while (1) {
         MTX_ASSERT_DEBUG(x >= 0 && x < deg);
         if (done[x])
            break;      // orbit complete
         done[x] = 1;
         x = p[x];
         ++orbitSize;
      }

      // Calculate the l.c.m of all orbit sizes
      order = lcm32u(order,orbitSize);
   }

   sysFree(done);

   return order;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
