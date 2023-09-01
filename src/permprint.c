////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Print a permutation
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

static int width(uint32_t x)
{
   if (x < 10) return 1;
   if (x < 100) return 2;
   if (x < 1000) return 3;
   if (x < 10000) return 4;
   if (x < 100000) return 5;
   if (x < 1000000) return 6;
   if (x < 10000000) return 7;
   if (x < 100000000) return 8;
   if (x < 1000000000) return 9;
   return 10;
}



/// @addtogroup perm
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Print a permutation
/// This function prints a permutation on the standard output using
/// cycle notation. If @em name is not 0, the name followed by an
/// equal sign is printed before the permutation. For example, the
/// statement <tt>permPrint("Perm",P);</tt> could produce the following output:
/// <pre>
/// Perm=(1 9)(2 3 6)(4 5 7)
/// </pre>
/// Fixed points are always suppressed in the output.
/// @param name Name to print before the permutation or 0.
/// @param perm Pointer to the permutation.
/// @return The function returns 0 on success and -1 on error.

void permPrint(const char *name, const Perm_t *perm)
{

   // check arguments
   permValidate(MTX_HERE, perm);

   // print the name
   if (name != NULL) {
      printf("%s=",name);
   }

   // print the permutation
   const uint32_t deg = perm->Degree;
   const uint32_t* p = perm->Data;

   uint8_t* done = NALLOC(uint8_t, deg);
   uint8_t* seed = done - 1;
   uint8_t* seedEnd = done + deg;

   // Run through all orbits. 
   int empty = 1;
   int count = 0;
   while (1) {
      // Find start of next orbit.
      while (*++seed != 0 && seed != seedEnd);
      if (seed == seedEnd)
         break;             // Done!
      uint32_t x = (uint32_t)(seed - done);

      // Suppress fixed points (GAP does not like them)
      if (p[x] == x) {
         done[x] = 1;
         continue;
      }

      empty = 0;
      int first = 1;
      while (1) {
         MTX_ASSERT_DEBUG(x >= 0 && x < deg);
         if (done[x])
            break;      // orbit complete
         done[x] = 1;

         if (first) {
            first = 0;
            if ((count += width(x) + 1) > 77) {
               printf("\n    (%lu",(unsigned long)x);
               count = 5 + width(x);
            } else {
               printf("(%lu",(unsigned long)x);
            }
         } else {
            if ((count += width(x)+1) > 77) {
               printf(",\n    %lu",(unsigned long)x);
               count = 4 + width(x);
            } else {
               printf(",%lu",(unsigned long)x);
            }
         }
         x = p[x];
      }
      printf(")");
      ++count;
   }

   sysFree(done);
   if (empty) {
      printf("()");
   }
   if (name != NULL) {
      printf("\n");
   }
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
