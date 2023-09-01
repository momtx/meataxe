////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Polynomial division
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup poly
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Polynomial division.
/// This function performs a polynomial division. Given two polynomials a and
/// b≠0 over the same field, %polDivMod() finds two  polynomials q and r such
/// that a=qb+r, and deg(r)<deg(b).
///
/// The quotient q is returned as the function result. This is a newly
/// allocated polynomial. The caller is responsible for deleting the quotient
/// when it no longer needed.
///
/// The remainder r, is stored in @em a and replaces the original value. If you
/// need to preserve the value of @em a you must make a copy using polDup() before
/// calling %polDivMod(). @em b is not changed.
/// @see polMod()
/// @param a First polynomial (numerator) on call, remainder on return.
/// @param b Second polynomial (denominator).
/// @return The quotient or 0 on error.

Poly_t *polDivMod(Poly_t *a, const Poly_t *b)
{
   Poly_t *q;

   // check arguments
   polValidate(MTX_HERE, a);
   polValidate(MTX_HERE, b);
   if (a->Field != b->Field) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
      return NULL;
   }
   ffSetField(a->Field);
   if (b->Degree <= -1) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_DIV0);
      return NULL;
   }
   if (a->Degree < b->Degree) {
      q = polAlloc(a->Field,-1);        // trivial case: Quotient = 0
   } else {
      FEL lead = b->Data[b->Degree];
      int i, k;

      if (lead == FF_ZERO) {
         mtxAbort(MTX_HERE,"%s",MTX_ERR_DIV0);
         return NULL;
      }
      q = polAlloc(ffOrder,a->Degree - b->Degree);
      if (q == NULL) {
         mtxAbort(MTX_HERE,"Cannot allocate result");
         return NULL;
      }
      for (i = a->Degree; i >= b->Degree; --i) {
         FEL qq = ffNeg(ffDiv(a->Data[i],lead));
         for (k = 0; k <= b->Degree; ++k) {
            a->Data[i - k] = ffAdd(a->Data[i - k],
                                   ffMul(qq,b->Data[b->Degree - k]));
         }
         q->Data[i - b->Degree] = ffNeg(qq);
      }
      Pol_Normalize(a);
   }
   return q;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Polynomial division.
/// This function replaces @em a with the remainder of the division of @em a by @em b.
/// @see polDivMod()
/// @param a First polynomial (numerator) on call, remainder on return.
/// @param b Second polynomial (denominator).
/// @return @em a or NULL on error.

Poly_t *polMod(Poly_t *a, const Poly_t *b)
{
   // check arguments
   polValidate(MTX_HERE,a);
   polValidate(MTX_HERE,b);
   if (a->Field != b->Field) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
      return NULL;
   }

   ffSetField(a->Field);
   if (b->Degree <= -1) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_DIV0);
      return NULL;
   }
   if (a->Degree >= b->Degree) {
      FEL lead = b->Data[b->Degree];
      int i, k;

      if (lead == FF_ZERO) {
         mtxAbort(MTX_HERE,"%s",MTX_ERR_DIV0);
         return NULL;
      }
      for (i = a->Degree; i >= b->Degree; --i) {
         FEL qq = ffNeg(ffDiv(a->Data[i],lead));
         for (k = 0; k <= b->Degree; ++k) {
            a->Data[i - k] = ffAdd(a->Data[i - k], ffMul(qq,b->Data[b->Degree - k]));
         }
	 MTX_ASSERT(a->Data[i] == FF_ZERO);
      }
      Pol_Normalize(a);
   }
   return a;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
