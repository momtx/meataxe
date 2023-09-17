////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Add two polynomials
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Extend buffer.

static void resize(Poly_t *p, int newdeg)
{
   int i;
   FEL *x;

   if (p->Degree < newdeg) {
      if (p->BufSize < newdeg + 1) {    // Allocate new buffer
         x = NREALLOC(p->Data,FEL,newdeg + 1);
         p->Data = x;
         p->BufSize = newdeg + 1;
      }
      for (i = p->Degree + 1; i <= newdeg; ++i) {
         p->Data[i] = FF_ZERO;
      }
      p->Degree = newdeg;
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// @addtogroup poly
/// @{
/// Add polynomials.
/// This function adds @em src to @em dest.
/// The polynomials must be over the same field.
/// @param dest First polynomial, overwritten with the result.
/// @param src Second polynomial.
/// @return @em dest, or 0 on error.

Poly_t *polAdd(Poly_t *dest, const Poly_t *src)
{
   FEL *s, *d;
   int i;

   polValidate(MTX_HERE, src);
   polValidate(MTX_HERE, dest);
   if (dest->Field != src->Field)
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
   if ((i = src->Degree) == -1) {
      return dest;      // src = 0 

   }
   ffSetField(src->Field);
   resize(dest,i);
   s = src->Data;
   d = dest->Data;
   for (; i >= 0; --i) {
      *d = ffAdd(*d,*s++);
      ++d;
   }
   Pol_Normalize(dest);
   return dest;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
