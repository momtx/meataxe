////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Print a factored polynomial
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup poly
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Multiply factored polynomials.
/// Multiplies @em dest by @em src. The previous content of @em dest is lost.
/// @see fpMulP()
/// @param dest Factored polynomial to modify.
/// @param src Factored polynomial.
/// @return The function returns |dest| or |NULL| on error.

FPoly_t *fpMul(FPoly_t *dest, const FPoly_t *src)
{
   int i;

   /* Check the arguments
      ------------------- */
   fpValidate(MTX_HERE, src);
   fpValidate(MTX_HERE, dest);

   for (i = 0; i < src->NFactors; ++i) {
      if (fpMulP(dest,src->Factor[i],src->Mult[i]) == NULL) {
         mtxAbort(MTX_HERE,"Cannot multiply");
         return NULL;
      }
   }
   return dest;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
