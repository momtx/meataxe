////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Print a factored polynomial.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

/*MTX_DEFINE_FILE_INFO*/

/// @addtogroup poly
/// @{

/// Print a factored polynomial.
/// This function prints a factored polynomial to the standard output.
/// If @em name is not 0, "name=" is printed before the polynomial.
/// @param name Name of the polynomial or 0.
/// @param p Pointer to the factored polynomial.
/// @return 0 on success, -1 on error.

int fpPrint(const char *name, const FPoly_t *p)
{
   int i;
   fpValidate(MTX_HERE, p);
   if (name != NULL) {
      printf("%s =",name);
   }
   for (i = 0; i < p->NFactors; ++i) {
      int e = p->Mult[i];
      if (i > 0) { printf("    * "); }
      printf("(");
      polPrint(NULL,p->Factor[i]);
      if (e > 1) {
         printf(")^%d\n",e);
      } else {
         printf(")\n");
      }
   }
   return 0;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
