////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Duplicate a factored polynomial
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

/// @addtogroup poly
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Creates a copy of a factored polynomial.

FPoly_t *fpDup(const FPoly_t *src)
{
   fpValidate(MTX_HERE, src);

   // Copy the factors
   Poly_t **new_factor = NALLOC(Poly_t *,src->NFactors);
   int* new_mult = NALLOC(int,src->NFactors);
   for (int i = 0; i < src->NFactors; ++i) {
      new_mult[i] = src->Mult[i];
      new_factor[i] = polDup(src->Factor[i]);
   }

   // Create a new factored polynomial
   FPoly_t *x = fpAlloc();
   sysFree(x->Factor);
   sysFree(x->Mult);
   x->Factor = new_factor;
   x->Mult = new_mult;
   x->BufSize = x->NFactors = src->NFactors;
   return x;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
