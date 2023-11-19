////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Print a factored polynomial.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

/// @addtogroup poly
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Format a factored polynomial

void fpFormat(StrBuffer* sb, const FPoly_t *p)
{
   fpValidate(MTX_HERE, p);
   if (p->NFactors == 0) {
      sbAppend(sb, "1");
      return;
   }

   for (int i = 0; i < p->NFactors; ++i) {
      int e = p->Mult[i];
      if (i > 0) { sbAppend(sb, " * "); }
      sbAppend(sb, "(");
      polFormat(sb, p->Factor[i]);
      if (e > 1) {
         sbPrintf(sb, ")^%d",e);
      } else {
         sbAppend(sb,")");
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////


/// Print a factored polynomial.
/// This function prints a factored polynomial to the standard output.
/// If @em name is not 0, "name=" is printed before the polynomial.
/// @param name Name of the polynomial or 0.
/// @param p Pointer to the factored polynomial.

void fpPrint(const char *name, const FPoly_t *p)
{
   fpValidate(MTX_HERE, p);
   StrBuffer* sb = sbAlloc(100);
   if (name != NULL) {
      sbPrintf(sb, "%s =",name);
   }
   fpFormat(sb, p);
   fputs(sbData(sb), stdout);
   sbFree(sb);
   if (name != NULL) {
      printf("\n");
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

char* fpToEphemeralString(const FPoly_t *p)
{
   StrBuffer* sb = sbAlloc(100);
   fpFormat(sb, p);
   return sbToEphemeralString(sb);
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
