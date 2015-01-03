////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Duplicate a factored polynomial
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
/// Duplicate a factored polynomial.
/// This function creates a copy of a factored polynomial.
/// @param src Pointer to a factored polynomial.
/// @return A pointer to a copy of @em src, or 0 on error.

FPoly_t *FpDup(const FPoly_t *src)
{
   FPoly_t *x;
   Poly_t **new_factor;
   int *new_mult;
   int i;

   /* Check the argument
      ------------------ */
   if (!FpIsValid(src)) {
      MTX_ERROR1("%E",MTX_ERR_BADARG);
      return NULL;
   }

   /* Copy the factors
      ---------------- */
   new_factor = NALLOC(Poly_t *,src->NFactors);
   if (new_factor == NULL) {
      MTX_ERROR1("%E",MTX_ERR_NOMEM);
      return NULL;
   }
   new_mult = NALLOC(int,src->NFactors);
   if (new_mult == NULL) {
      MTX_ERROR1("%E",MTX_ERR_NOMEM);
      return NULL;
   }
   for (i = 0; i < src->NFactors; ++i) {
      new_mult[i] = src->Mult[i];
      new_factor[i] = PolDup(src->Factor[i]);
      if (new_factor[i] == NULL) {
         while (--i >= 0) {
            PolFree(new_factor[i]);
         }
         SysFree(new_factor);
         SysFree(new_mult);
         MTX_ERROR1("Cannot duplicate polynomial",MTX_ERR_NOMEM);
         return NULL;
      }
   }

   /* Create a new factored polynomial
      ------------------------------- */
   x = FpAlloc();
   if (x == NULL) {
      MTX_ERROR("Cannot create copy");
      return NULL;
   }
   SysFree(x->Factor);
   SysFree(x->Mult);
   x->Factor = new_factor;
   x->Mult = new_mult;
   x->BufSize = x->NFactors = src->NFactors;
   return x;
}


/// @}