////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Polynomial division
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

/// @addtogroup poly
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Polynomial division.
/// This function performs a polynomial division. Given two polynomials a and
/// b≠0 over the same field, %PolDivMod() finds two  polynomials q and r such
/// that a=qb+r, and deg(r)<deg(b).
///
/// The quotient q is returned as the function result. This is a newly
/// allocated polynomial. The caller is responsible for deleting the quotient
/// when it no longer needed.
///
/// The remainder r, is stored in @em a and replaces the original value. If you
/// need to preserve the value of @em a you must make a copy using PolDup() before
/// calling %PolDivMod(). @em b is not changed.
/// @see PolMod()
/// @param a First polynomial (numerator) on call, remainder on return.
/// @param b Second polynomial (denominator).
/// @return The quotient or 0 on error.

Poly_t *PolDivMod(Poly_t *a, const Poly_t *b)
{
   Poly_t *q;

   // check arguments
   if (!PolIsValid(a) || !PolIsValid(b)) {
      return NULL;
   }
   if (a->Field != b->Field) {
      MTX_ERROR1("%E",MTX_ERR_INCOMPAT);
      return NULL;
   }
   FfSetField(a->Field);
   if (b->Degree <= -1) {
      MTX_ERROR1("%E",MTX_ERR_DIV0);
      return NULL;
   }
   if (a->Degree < b->Degree) {
      q = PolAlloc(a->Field,-1);        // trivial case: Quotient = 0
   } else {
      FEL lead = b->Data[b->Degree];
      int i, k;

      if (lead == FF_ZERO) {
         MTX_ERROR1("%E",MTX_ERR_DIV0);
         return NULL;
      }
      q = PolAlloc(FfOrder,a->Degree - b->Degree);
      if (q == NULL) {
         MTX_ERROR("Cannot allocate result");
         return NULL;
      }
      for (i = a->Degree; i >= b->Degree; --i) {
         FEL qq = FfNeg(FfDiv(a->Data[i],lead));
         for (k = 0; k <= b->Degree; ++k) {
            a->Data[i - k] = FfAdd(a->Data[i - k],
                                   FfMul(qq,b->Data[b->Degree - k]));
         }
         q->Data[i - b->Degree] = FfNeg(qq);
      }
      Pol_Normalize(a);
   }
   return q;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Polynomial division.
/// This function replaces @em a with the remainder of the division of @em a by @em b.
/// @see PolDivMod()
/// @param a First polynomial (numerator) on call, remainder on return.
/// @param b Second polynomial (denominator).
/// @return @em a or 0 on error.

Poly_t *PolMod(Poly_t *a, const Poly_t *b)
{
   // check arguments
   if (!PolIsValid(a) || !PolIsValid(b)) {
      return NULL;
   }
   if (a->Field != b->Field) {
      MTX_ERROR1("%E",MTX_ERR_INCOMPAT);
      return NULL;
   }

   FfSetField(a->Field);
   if (b->Degree <= -1) {
      MTX_ERROR1("%E",MTX_ERR_DIV0);
      return NULL;
   }
   if (a->Degree >= b->Degree) {
      FEL lead = b->Data[b->Degree];
      int i, k;

      if (lead == FF_ZERO) {
         MTX_ERROR1("%E",MTX_ERR_DIV0);
         return NULL;
      }
      for (i = a->Degree; i >= b->Degree; --i) {
         FEL qq = FfNeg(FfDiv(a->Data[i],lead));
         for (k = 0; k <= b->Degree; ++k) {
            a->Data[i - k] = FfAdd(a->Data[i - k], FfMul(qq,b->Data[b->Degree - k]));
         }
      }
      Pol_Normalize(a);
   }
   return a;
}


/// @}
