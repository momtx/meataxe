////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Basic polynomial functions.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

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
   if (p == NULL || p->typeId != MTX_TYPE_POLYNOMIAL) {
      return 0;
   }
   const int deg = p->degree;
   if (deg < -1 || p->field < 2 || p->data == NULL || p->bufSize < 0) {
      return 0;
   }
   if (deg >= 0 && p->data[deg] == FF_ZERO) {
      return 0;
   }
   return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Checks if the given polynomial is valid and aborts the program if the test fails.

void polValidate(const struct MtxSourceLocation* src, const Poly_t* pol)
{
   if (pol == NULL) {
      mtxAbort(src, "NULL polynomial");
   }
   if (pol->typeId != MTX_TYPE_POLYNOMIAL || pol->degree < -1 || pol->field < 2) {
      mtxAbort(src, "Invalid polynomial (typeId=0x%lx, field=%lu, deg=%ld)",
               (unsigned long) pol->typeId, (unsigned long) pol->field, (long) pol->degree);
   }
   if (pol->data == NULL || pol->bufSize < 0) {
      mtxAbort(src, "Invalid polynomial (data=0x%lx, size=%d)",
               (unsigned long)pol->data, pol->bufSize);
   }
   if (pol->degree >= 0 && pol->data[pol->degree] == FF_ZERO) {
      mtxAbort(src, "Invalid polynomial (leading coefficient is zero)");
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Allocate a polynomial
/// This function creates the polynomial p(x)=x^n over the current field.
/// If n is negative, a zero polynomial (degree -1) is created. The return value is a pointer to a
/// newly allocated Poly_t structure. The caller is responsible for releasing memory by calling
/// @ref polFree when the polynomial is no longer needed.
/// @param field Field order.
/// @param degree Degree of the polynomial.
/// @return Pointer to a new Poly_t structure

Poly_t* polAlloc(uint32_t field, int32_t degree)
{
   if (degree < 0) {degree = -1;}
   MTX_ASSERT(degree < (int32_t) 2147483647L);

   ffSetField(field);
   Poly_t* x = ALLOC(Poly_t);
   x->typeId = MTX_TYPE_POLYNOMIAL;
   x->field = field;
   x->degree = degree;
   x->bufSize = degree + 1;
   x->data = NALLOC(FEL, x->bufSize);
   for (int32_t i = 0; i < degree; ++i) {
      x->data[i] = FF_ZERO;
   }
   if (degree >= 0) {
      x->data[degree] = FF_ONE;
   }
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
   sysFree(x->data);
   memset(x,0,sizeof(Poly_t));
   sysFree(x);
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Normalize a polynomial.
/// This function makes sure that the leading coefficient of a polynomial is non-zero.

void Pol_Normalize(Poly_t *p)
{
   int i = p->degree;
   while (i >= 0 && p->data[i] == FF_ZERO) {
      --i;
   }
   p->degree = i;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
