////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Multiply two polynomials
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

/// @addtogroup poly
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Multiply polynomials.
/// This function multiplies @em dest by @em src and returns @em dest.
/// The polynomials must be over the same field.
/// @param dest Pointer to the first polynomial.
/// @param src Pointer to the second polynomial.
/// @return @em dest, or 0 on error.

Poly_t *PolMul(Poly_t *dest, const Poly_t *src)
{
   FEL *x, *y, *d, *s;
   int di, si;
   size_t xdeg;

   // check arguments
   if (!PolIsValid(src) || !PolIsValid(dest)) {
      return NULL;
   }
   if (dest->Field != src->Field) {
      MTX_ERROR1("%E",MTX_ERR_INCOMPAT);
      return NULL;
   }

   // handle special cases: dest = 0, src = 0
   if (dest->Degree == -1) {
      return dest;
   }
   if (src->Degree == -1) {
      dest->Degree = -1;
      return dest;
   }

   d = dest->Data;
   s = src->Data;
   xdeg = src->Degree + dest->Degree;
   FfSetField(src->Field);

   // allocate result buffer
   x = NALLOC(FEL,xdeg + 1);
   if (x == NULL) {
      MTX_ERROR("Cannot allocate result");
      return NULL;
   }
   for (di = xdeg, y = x; di >= 0; --di) {
      *y++ = FF_ZERO;
   }

   // multiply
   for (di = 0; di <= dest->Degree; ++di) {
      for (si = 0; si <= src->Degree; ++si) {
         x[si + di] = FfAdd(x[si + di],FfMul(s[si],d[di]));
      }
   }

   // overwrite <dest> with the result
   SysFree(dest->Data);
   dest->Data = x;
   dest->Degree = xdeg;
   dest->BufSize = xdeg + 1;
   return dest;
}


/// @}
