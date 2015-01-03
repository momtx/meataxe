////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Multiply factored polynomial
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

/// @addtogroup poly
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Multiply with an irreducible polynomial.
/// This function multiplies a factored polynomial with the power of an
/// an irreducible factor. It is not checked that @em src is irreducible.
/// @see FpMul()
/// @param dest Factored polynomial to modify.
/// @param src Irreducible polynomial.
/// @param pwr Power of the irreducible polynomial.
/// @return The function returns @em dest or 0 on error.

FPoly_t *FpMulP(FPoly_t *dest, const Poly_t *src, int pwr)
{
   int i;
   int cmp = 0;

   /* Check the arguments
      ------------------- */
   if (!PolIsValid(src) || !FpIsValid(dest)) {
      return NULL;
   }
   if (pwr <= 0) {
      MTX_ERROR2("pwr=%d: %E",pwr,MTX_ERR_BADARG);
      return NULL;
   }

   /* Find the insert position
      ------------------------ */
   for (i = 0;
        i < dest->NFactors && (cmp = PolCompare(dest->Factor[i],src)) < 0;
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

         if ((e == NULL) || (x == NULL)) {
            MTX_ERROR("Cannot grow: %S");
            return NULL;
         }
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
      dest->Factor[i] = PolDup(src);
      dest->Mult[i] = pwr;
      if (dest->Factor[i] == NULL) {
         MTX_ERROR("Cannot copy polynomial");
         return NULL;
      }
   } else {
      dest->Mult[i] += pwr;
   }
   return dest;
}


/// @}