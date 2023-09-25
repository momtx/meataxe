////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Add two polynomials
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

/// @addtogroup poly
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

static void grow(Poly_t *p, int newdeg)
{
   if (p->degree < newdeg) {
      if (p->bufSize < newdeg + 1) {
         p->bufSize = newdeg + 1;
         p->data = NREALLOC(p->data, FEL, p->bufSize);
      }
      for (int32_t i = p->degree + 1; i <= newdeg; ++i) {
         p->data[i] = FF_ZERO;
      }
      p->degree = newdeg;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Add polynomials.
/// This function adds @em src to @em dest.
/// The polynomials must be over the same field.
/// @param dest First polynomial, overwritten with the result.
/// @param src Second polynomial.
/// @return @em dest, or 0 on error.

Poly_t *polAdd(Poly_t *dest, const Poly_t *src)
{
   polValidate(MTX_HERE, src);
   polValidate(MTX_HERE, dest);
   if (dest->field != src->field)
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
   if (src->degree == -1) {
      return dest;      // src = 0 
   }
   ffSetField(src->field);
   grow(dest,src->degree);
   FEL* s = src->data;
   FEL* d = dest->data;
   for (int32_t i = src->degree; i >= 0; --i) {
      *d = ffAdd(*d,*s++);
      ++d;
   }
   Pol_Normalize(dest);
   return dest;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
