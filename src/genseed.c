////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Seed vector generator
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <limits.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Local data


////////////////////////////////////////////////////////////////////////////////////////////////////
/// @defgroup pseed The seed vector generator
///
/// The seed vector generator is used to walk through the one-dimensional subspaces of a given
/// vector space V, the @em "seed space". For each one-dimensional subspace U≤V the generator
/// produces a representant u∊U. These vectors are called @em "seed vectors".
///
/// Once a basis b<sub>1</sub>,...b<sub>n</sub> for the seed space is fixed, each vector
/// v=λ<sub>1</sub>b<sub>1</sub>+...+λ<sub>n</sub>b<sub>n</sub> can be identified by a natural
/// number by mapping coefficients to natural numbers in the usual way (see ffToInt() for details)
/// and calculating λ<sub>1</sub>q<sup>0</sup>+...λ<sub>n</sub>q<sup>n-1</sup>.
/// Seed vectors are those vectors where the leading coefficient is 1.

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @ingroup pseed
/// Generates the next seed vector.
/// @a lastno is the number of the previous seed vector and can be zero to calculate the first seed
/// vector.
/// The seed vector is stored in @a vec, and its number is returned.
/// If the rows of @a basis are not linearly independent, there will be redundant seed
/// vectors, but no error occurs.
///
/// Typically, MakeSeedVector() is called repeatedly until the whole
/// seed space is exhausted. Here is s short example:
/// @code
/// long v;
/// PTR vec = ffAlloc(1);
/// v = MakeSeedVector(basis,0,ptr);
/// while (v > 0)
/// {
///     ...
///     v = MakeSeedVector(basis,v,ptr)
/// }
/// @endcode
///
/// @param basis Basis of the seed space. Need not be in echelon form.
/// @param lastno Previous seed vector (0 to start from beginning).
/// @param vec Buffer for the seed vector.
/// @return Seed vector number, or -1 on error.

long MakeSeedVector(const Matrix_t *basis, long lastno, PTR vec)
{
   long nextno, x, i;
   int j;

   matValidate(MTX_HERE, basis);
   if ((vec == NULL) || (lastno < 0)) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_BADARG);
      return -1;
   }

   // Find the next seed vector number
   nextno = lastno + 1;
   for (x = 1; (i = nextno / x) >= ffOrder; x *= ffOrder) {
   }
   if (i != 1) {
      nextno = x * ffOrder;
   }

   // Make the seed vector
   ffSetField(basis->Field);
   ffMulRow(vec,FF_ZERO, basis->Noc);
   for (j = 0, x = nextno; x != 0 && j < basis->Nor; ++j, x /= ffOrder) {
      FEL co = ffFromInt(x % ffOrder);
      if (co != FF_ZERO) {
         ffAddMulRow(vec,matGetPtr(basis,j),co, basis->Noc);
      }
   }

   // Check for overflow
   if (x != 0) {
      return -1;
   }

   return nextno;
}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
