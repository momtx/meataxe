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
   if (a->field != b->field) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
      return NULL;
   }
   ffSetField(a->field);
   if (b->degree <= -1) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_DIV0);
      return NULL;
   }
   if (a->degree < b->degree) {
      q = polAlloc(a->field,-1);        // trivial case: Quotient = 0
   } else {
      FEL lead = b->data[b->degree];
      int i, k;

      if (lead == FF_ZERO) {
         mtxAbort(MTX_HERE,"%s",MTX_ERR_DIV0);
         return NULL;
      }
      q = polAlloc(ffOrder,a->degree - b->degree);
      if (q == NULL) {
         mtxAbort(MTX_HERE,"Cannot allocate result");
         return NULL;
      }
      for (i = a->degree; i >= b->degree; --i) {
         FEL qq = ffNeg(ffDiv(a->data[i],lead));
         for (k = 0; k <= b->degree; ++k) {
            a->data[i - k] = ffAdd(a->data[i - k],
                                   ffMul(qq,b->data[b->degree - k]));
         }
         q->data[i - b->degree] = ffNeg(qq);
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
   if (a->field != b->field) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
      return NULL;
   }

   ffSetField(a->field);
   if (b->degree <= -1) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_DIV0);
      return NULL;
   }
   if (a->degree >= b->degree) {
      FEL lead = b->data[b->degree];
      int i, k;

      if (lead == FF_ZERO) {
         mtxAbort(MTX_HERE,"%s",MTX_ERR_DIV0);
         return NULL;
      }
      for (i = a->degree; i >= b->degree; --i) {
         FEL qq = ffNeg(ffDiv(a->data[i],lead));
         for (k = 0; k <= b->degree; ++k) {
            a->data[i - k] = ffAdd(a->data[i - k], ffMul(qq,b->data[b->degree - k]));
         }
	 MTX_ASSERT(a->data[i] == FF_ZERO);
      }
      Pol_Normalize(a);
   }
   return a;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
