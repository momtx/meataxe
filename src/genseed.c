////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Seed vector generator
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <inttypes.h>
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
/// number by mapping the coefficients to natural numbers in the usual way (see ffToInt() for
/// details) and treating them as as digits in a base q representation of the vector number.
/// Seed vectors are those vectors where the leading digit (after erasing leading zeroes) is 1.

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @ingroup pseed
/// Calculates a seed vector given its number. See also @ref svgMakeNext.
/// The function fails and aborts the program if @a number is not a valid seed vector number.
/// In particular, passing @a number = 0 will always fail.
///
/// @param vec Buffer for the seed vector.
/// @param number The seed vector number.
/// @param basis Basis of the seed space. Need not be in echelon form.

void svgMake(PTR vec, uint32_t number, const Matrix_t *basis)
{
   matValidate(MTX_HERE, basis);
   if (vec == NULL || number < 0) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_BADARG);
   }
   ffSetField(basis->field);
   ffMulRow(vec,FF_ZERO, basis->noc);
   uint32_t x = number;
   for (uint32_t j = 0; x != 0 && j < basis->nor; ++j) {
      FEL co = ffFromInt(x % ffOrder);
      if (co != FF_ZERO) {
         ffAddMulRow(vec,matGetPtr(basis,j),co, basis->noc);
      }
      x /= ffOrder;
   }

   if (x != 0)
      mtxAbort(MTX_HERE,"Bad seed vector number %"PRIu32, number);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @ingroup pseed
/// Generates the next seed vector. Returns 0 on success and -1 if there are no more seed vectors.
///
/// When the first seed vector is calculated (i.e., if @a number points to a zero value) the
/// function verifies that all seed vector numbers can be represented as 32 bit unsigned integers. 
/// If this is not the case, the function fails and aborts the program. Note that this check is
/// only performed for the first seed vector (with number 1).
/// In all other cases, the function will fail only if the next vector number would be greater than
/// 2<sup>32</sup>-1
///
/// @a basis is the basis for the seed space. It is not checked whether the bases vectors are
/// linearly independent. If they are not, there will be redundant seed vectors but not no error
/// occurs.
///
/// @a vec is a row buffer which is overwritten with the generated seed vector.
/// @a vec may be NULL if the seed vector is not required.
///
/// @a number points to a variable containing the previous seed vector number (or 0 if no seed
/// vector has been calculated yet. After successful return, the variable is updated and contains
/// the number of the generated seed vector. If svgNext() fails, the variable is not changed
/// (however, @a vec may have been be modified).
///
/// Usage example:
/// @code
/// Matrix_t* basis = ...;
/// long vecno = 0;
/// PTR vec = ffAlloc(1, basis->noc);
/// while (svgNext(basis, &vecno, vec) == 0) {
/// {
///     spinUpWithSeed(vec, ...);
/// }
/// @endcode

int svgMakeNext(PTR vec, uint32_t* number, const Matrix_t* basis)
{
   matValidate(MTX_HERE, basis);
   if (number == NULL) {
      mtxAbort(MTX_HERE, "%s", MTX_ERR_BADARG);
   }
   if (*number == 0) {
      uint64_t x = 1;
      for (uint32_t i = 1; i < basis->nor; ++i) {
         x = (x + (ffOrder - 1)) * ffOrder;
         if (x >= 0x100000000ULL) {
            mtxAbort(
               MTX_HERE,
               "Seed space to large (q=%"PRIu32" nor=%"PRIu32")", ffOrder, basis->nor);
         }
      }
   }

   // Find the next seed vector number (leading digit in base q must be 1)
   uint32_t nextNumber = *number + 1;
   {
      uint32_t x = 1;
      while (x * ffOrder < nextNumber) {
         x *= ffOrder;
      }
      if (nextNumber >= 2 * x) { nextNumber = x * ffOrder; }
   }

   // Make the seed vector (if requested)
   ffSetField(basis->field);
   uint32_t x = nextNumber;
   if (vec != NULL) {
      ffMulRow(vec, FF_ZERO, basis->noc);
      for (unsigned j = 0; x != 0 && j < basis->nor; ++j) {
         FEL co = ffFromInt(x % ffOrder);
         if (co != FF_ZERO) {
            ffAddMulRow(vec, matGetPtr(basis, j), co, basis->noc);
         }
         x /= ffOrder;
      }
   } else {
      for (unsigned j = 0; x != 0 && j < basis->nor; ++j)
         x /= ffOrder;
   }

   // Check if the vector number is valid
   if (x != 0) {
      return -1;
   } else {
      *number = nextNumber;
      return 0;
   }
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
