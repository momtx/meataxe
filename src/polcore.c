////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Basic polynomial functions.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


#define POLY_MAGIC 0x355A3207

/// @defgroup poly Polynomials
/// @details
/// The MeatAxe can work with polynomials over a finite field. A polynomial is represented
/// by a Poly_t structure. Each polynomial carries the field order, i.e., you can work
/// with polynomials over different fields on one program. However, this feature is
/// currently of little use since all standard operations only work on polynomials over the
/// same field, and there is no easy way to identify polynomials over a field and its subfields.
///
/// There is a second representation of polynomials as product of factors, see FPoly_t.

/// @addtogroup poly
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @class Poly_t
/// A polynomial over a finite field.
/// Internally, a polynomial of degree n is represented as an array of n+1 field
/// elements (@c Data field), where <tt>data[i]</tt> is the coefficient of x^i.
/// The leading coefficient is always non-zero on the MeatAxe API level (it can
/// temporarily be zero during calculations). The null polynomial has a degree of -1.

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Check a polynomial.
/// This function checks if the argument is a pointer to a valid polynomial. If the polynomial
/// the function returns 1. Otherwise, an error is signalled and, if the error handler does not
/// terminate the program, the function returns 0.
/// @param p The polynomial to check.
/// @return 1 if @em p points to a valid polynomial, 0 otherwise.

int polIsValid(const Poly_t *p)
{
   if (p == NULL || p->Magic != POLY_MAGIC) {
      return 0;
   }
   const int deg = p->Degree;
   if (deg < -1 || p->Field < 2 || p->Data == NULL || p->BufSize < 0) {
      return 0;
   }
   if (deg >= 0 && p->Data[deg] == FF_ZERO) {
      return 0;
   }
   return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Checks if the given polynimial is valid and aborts the program if the test fails.

void polValidate(const struct MtxSourceLocation* src, const Poly_t *pol)
{
   if (pol == NULL)
      mtxAbort(src,"NULL polynomial");
   if (pol->Magic != POLY_MAGIC || pol->Degree < -1 || pol->Field < 2) {
      mtxAbort(src,"Invalid polynomial (magic=0x%x, field=%d, deg=%d)",
            pol->Magic, pol->Field, pol->Degree);
   }
   if (pol->Data == NULL || pol->BufSize < 0) {
      mtxAbort(src,"Invalid polynomial (data=0x%lx, size=%d)",
            (unsigned long)pol->Data, pol->BufSize);
   }
   if (pol->Degree >= 0 && pol->Data[pol->Degree] == FF_ZERO) {
      mtxAbort(src,"Invalid polynomial (leading coefficient is zero)");
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Allocate a polynomial
/// This function creates the polynomial p(x)=x^n over the current field.
/// If n is negative, a zero polynomial is created. The return value is a
/// pointer to a newly allocated Poly_t structure. The caller is responsible
/// for releasing memory by calling polFree() when the polynomial is no
/// longer needed.
/// @param fl Field order.
/// @param n Degree of the polynomial.
/// @return Pointer to a new Poly_t structure

Poly_t *polAlloc(int fl, int n)
{
   Poly_t *x;
   int i, s;

   if (n < 0) {
      n = -1;
   }
   if ((s = n + 1) < 1) {
      s = 1;
   }

   ffSetField(fl);
   if ((x = ALLOC(Poly_t)) == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate polynomial");
   }
   x->Magic = POLY_MAGIC;
   x->Field = fl;
   x->Degree = n;
   x->BufSize = s;
   if ((x->Data = NALLOC(FEL,s)) == NULL) {
      sysFree(x);
      mtxAbort(MTX_HERE,"Cannot allocate polynomial data");
   }
   for (i = 0; i < (int) s - 1; ++i) {
      x->Data[i] = FF_ZERO;
   }
   x->Data[s - 1] = FF_ONE;
   return x;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Free a polynomial
/// This function frees a polynomial data structure and cleans up all internal data.
/// @param x Pointer to the polynomial.
/// @return $0$ on success, $-1$ on error.

int polFree(Poly_t *x)
{
   polValidate(MTX_HERE, x);
   sysFree(x->Data);
   memset(x,0,sizeof(Poly_t));
   sysFree(x);
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Normalize a polynomial.
/// This function makes sure that the leading coefficient of a polynomial is non-zero.

void Pol_Normalize(Poly_t *p)
{
   int i = p->Degree;
   while (i >= 0 && p->Data[i] == FF_ZERO) {
      --i;
   }
   p->Degree = i;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
