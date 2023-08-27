////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Compare polynomials
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup poly
/// @{

// Gives identical results for ZZZ=0 and ZZZ=1.
static inline int ffCompare(FEL a, FEL b)
{
#if MTX_ZZZ == 0
   return a > b ? 1 : ((a == b) ?  0 : -1);
#elif MTX_ZZZ == 1
   int d = ffToInt(a) - ffToInt(b);
   return d > 0 ? 1 : (d == 0 ?  0 : -1);
#else
   #error
#endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Compare Two poynomials.
/// This function compares two polynomials and returns 0 if the polynomials are equal,
/// -1 if a<b, or 1 if a>b.  The ordering of polynomials is defined as follows.
/// - If a and b are over different fields, the polynomials over
///   the larger field is greater.
/// - Otherwise, if they have different degrees, the polynomial
///   with the higher degree is greater.
/// - If both field and degree are equal, the result of the comparison is 0
///   if the polynomials are equal. Otherwise it is unspecified if the return value
///   is +1 or -1.
/// @param a First polynomial.
/// @param b Second polynomial.
/// @return 1 if a>b, -1 if a<b, 0 if a=b, or -2 on error.

int polCompare(const Poly_t *a, const Poly_t *b)
{
   int i;

   // check arguments
   polValidate(MTX_HERE, a);
   polValidate(MTX_HERE, b);

   // compare fields
   if (a->Field > b->Field)
      return 1;
   if (a->Field < b->Field)
      return -1;

   // compare degrees
   if (a->Degree > b->Degree) {
      return 1;
   }
   if (a->Degree < b->Degree) {
      return -1;
   }

   // compare coefficients
   for (i = a->Degree; i >= 0; --i) {
       int cmp = ffCompare(a->Data[i], b->Data[i]);
       if (cmp != 0)
	   return cmp;
   }

   // the polynomials are equal
   return 0;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
