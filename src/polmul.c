////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Multiply two polynomials
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup poly
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Multiply polynomials.
/// This function multiplies @em dest by @em src and returns @em dest.
/// The polynomials must be over the same field.
/// @param dest Pointer to the first polynomial.
/// @param src Pointer to the second polynomial.
/// @return @em dest, or 0 on error.

Poly_t *polMul(Poly_t *dest, const Poly_t *src)
{
   FEL *x, *y, *d, *s;
   int di, si;
   size_t xdeg;

   // check arguments
   polValidate(MTX_HERE, src);
   polValidate(MTX_HERE, dest);
   if (dest->field != src->field) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
      return NULL;
   }

   // handle special cases: dest = 0, src = 0
   if (dest->degree == -1) {
      return dest;
   }
   if (src->degree == -1) {
      dest->degree = -1;
      return dest;
   }

   d = dest->data;
   s = src->data;
   xdeg = src->degree + dest->degree;
   ffSetField(src->field);

   // allocate result buffer
   x = NALLOC(FEL,xdeg + 1);
   if (x == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate result");
      return NULL;
   }
   for (di = xdeg, y = x; di >= 0; --di) {
      *y++ = FF_ZERO;
   }

   // multiply
   for (di = 0; di <= dest->degree; ++di) {
      for (si = 0; si <= src->degree; ++si) {
         x[si + di] = ffAdd(x[si + di],ffMul(s[si],d[di]));
      }
   }

   // overwrite <dest> with the result
   sysFree(dest->data);
   dest->data = x;
   dest->degree = xdeg;
   dest->bufSize = xdeg + 1;
   return dest;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
