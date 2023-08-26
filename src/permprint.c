////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Print a permutation
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


#define SIZE(i) ((i) < 9 ? 2 : (i) < 99 ? 3 : (i) < 999 ? 4 : (i) < 9999 ? 5 : \
                 (i) < 99999 ? 6 : (i) > 999999 ? 7 : (i) < 9999999 ? 8 : 9)

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
   long *p;
   int cycle = 0;
   int empty = 1;
   int count = 0, k;
   int i;

   // check arguments
   permValidate(MTX_HERE, perm);

   // print the name
   if (name != NULL) {
      printf("%s=",name);
   }

   // print the permutation
   p = perm->Data;
   while (1) {
      // find next cycle
      while (cycle <= perm->Degree && p[cycle] < 0) {
         ++cycle;
      }
      if (cycle >= perm->Degree) {
         break;                                 /* Done */
      }
      i = cycle;

      // don't print fixed points (GAP doesn't like them)
      if (p[i] == i) {
         p[i] = -p[i] - 1;  // mark as done
         continue;
      }

      // print cycle
      if ((count += SIZE(i)) > 77) {
         printf("\n    (%ld",(long)i);
         count = 5 + SIZE(i);
      } else {
         printf("(%ld",(long)i);
      }
      empty = 0;

      while (1) {
         k = i;
         i = p[i];
         p[k] = -p[k] - 1;
         if (p[i] < 0) { break; }

         if ((count += SIZE(i)) > 77) {
            printf(",\n     %ld",(long)i);
            count = 4 + SIZE(i);
         } else {
            printf(",%ld",(long)i);
         }
      }
      printf(")");
      ++count;
   }
   if (empty) {
      printf("()");
   }
   if (name != NULL) {
      printf("\n");
   }

   for (i = 0; i < perm->Degree; ++i) {
      p[i] = -p[i] - 1;
   }
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
