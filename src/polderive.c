////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Polynomial derivation
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup poly
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Derive a polynomial.
/// This function derives a polynomial. Note that the derived polynomial is
/// stored in @em pol, replacing the original polynomial. The following piece of
/// code shows how to keep the original polynomial intact while calculating
/// the derivative:
/// @code
/// Poly_t *pol, *der;
/// ...
/// der = polDerive(polDup(pol));
/// @endcode
/// @param pol Pointer to the polynomial.
/// @return @em pol.

Poly_t *polDerive(Poly_t *pol)
{
   int i, maxdeg = -1;
   register FEL *buf;
   FEL f = FF_ZERO;

   // check argument
   polValidate(MTX_HERE, pol);

   buf = pol->Data;
   ffSetField(pol->Field);
   for (i = 0; i < pol->Degree; ++i) {
      f = ffAdd(f,FF_ONE);
      buf[i] = ffMul(buf[i + 1],f);
      if (buf[i] != FF_ZERO) {
         maxdeg = i;
      }
   }
   pol->Degree = maxdeg;
   return pol;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
