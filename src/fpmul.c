////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Multiply factored polynomial
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup poly
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Multiply with an irreducible polynomial.
/// This function multiplies a factored polynomial with the power of an
/// an irreducible factor. It is not checked that @em src is irreducible.
/// @see fpMul()
/// @param dest Factored polynomial to modify.
/// @param src Irreducible polynomial.
/// @param pwr Power of the irreducible polynomial.
/// @return The function returns @em dest or 0 on error.

FPoly_t *fpMulP(FPoly_t *dest, const Poly_t *src, int pwr)
{
   int i;
   int cmp = 0;

   /* Check the arguments
      ------------------- */
   polValidate(MTX_HERE, src);
   fpValidate(MTX_HERE, dest);
   if (pwr <= 0) {
      mtxAbort(MTX_HERE,"pwr=%d: %s",pwr,MTX_ERR_BADARG);
      return NULL;
   }

   /* Find the insert position
      ------------------------ */
   for (i = 0;
        i < dest->NFactors && (cmp = polCompare(dest->Factor[i],src)) < 0;
        ++i) {
   }

   /* Extend the buffer, if necessary
      ------------------------------- */
   if ((i >= dest->NFactors) || (cmp != 0)) {
      int k;
      if (dest->NFactors >= dest->BufSize) {
         int newsize = dest->BufSize + 5;
         Poly_t **x = NREALLOC(dest->Factor,Poly_t *,newsize);
         int *e = NREALLOC(dest->Mult,int,newsize);
         dest->Factor = x;
         dest->Mult = e;
         dest->BufSize = newsize;
      }

      /* Make room for the new factor
         ---------------------------- */
      for (k = dest->NFactors; k > i; --k) {
         dest->Factor[k] = dest->Factor[k - 1];
         dest->Mult[k] = dest->Mult[k - 1];
      }
      ++dest->NFactors;

      /* Insert new factor
         ----------------- */
      dest->Factor[i] = polDup(src);
      dest->Mult[i] = pwr;
      if (dest->Factor[i] == NULL) {
         mtxAbort(MTX_HERE,"Cannot copy polynomial");
         return NULL;
      }
   } else {
      dest->Mult[i] += pwr;
   }
   return dest;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
