////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Duplicate a factored polynomial
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup poly
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Duplicate a factored polynomial.
/// This function creates a copy of a factored polynomial.
/// @param src Pointer to a factored polynomial.
/// @return A pointer to a copy of @em src, or 0 on error.

FPoly_t *fpDup(const FPoly_t *src)
{
   FPoly_t *x;
   Poly_t **new_factor;
   int *new_mult;
   int i;

   /* Check the argument
      ------------------ */
   fpValidate(MTX_HERE, src);

   /* Copy the factors
      ---------------- */
   new_factor = NALLOC(Poly_t *,src->NFactors);
   if (new_factor == NULL) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_NOMEM);
      return NULL;
   }
   new_mult = NALLOC(int,src->NFactors);
   if (new_mult == NULL) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_NOMEM);
      return NULL;
   }
   for (i = 0; i < src->NFactors; ++i) {
      new_mult[i] = src->Mult[i];
      new_factor[i] = polDup(src->Factor[i]);
      if (new_factor[i] == NULL) {
         while (--i >= 0) {
            polFree(new_factor[i]);
         }
         sysFree(new_factor);
         sysFree(new_mult);
         mtxAbort(MTX_HERE,"Cannot duplicate polynomial",MTX_ERR_NOMEM);
         return NULL;
      }
   }

   /* Create a new factored polynomial
      ------------------------------- */
   x = fpAlloc();
   if (x == NULL) {
      mtxAbort(MTX_HERE,"Cannot create copy");
      return NULL;
   }
   sysFree(x->Factor);
   sysFree(x->Mult);
   x->Factor = new_factor;
   x->Mult = new_mult;
   x->BufSize = x->NFactors = src->NFactors;
   return x;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
