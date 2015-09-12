////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Seed vector generator
//
// (C) Copyright 1998-2014 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>
#include <limits.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Local data

MTX_DEFINE_FILE_INFO

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @defgroup pseed The seed vector generator
///
/// The seed vector generator is used to walk through the one-dimensional subspaces of a given
/// vector space V, the @em "seed space". For each one-dimensional subspace U≤V the generator
/// produces a representant u∊U. These vectors are called @em "seed vectors".
///
/// Once a basis b<sub>1</sub>,...b<sub>n</sub> for the seed space is fixed, each vector
/// v=λ<sub>1</sub>b<sub>1</sub>+...+λ<sub>n</sub>b<sub>n</sub> can be identified by a natural
/// number by mapping coefficients to natural numbers in the usual way (see FfToInt() for details)
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
/// PTR vec = FfAlloc(1);
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

   if (!MatIsValid(basis)) {
      return -1;
   }
   if ((vec == NULL) || (lastno < 0)) {
      MTX_ERROR1("%E",MTX_ERR_BADARG);
      return -1;
   }

   // Find the next seed vector number
   nextno = lastno + 1;
   for (x = 1; (i = nextno / x) >= FfOrder; x *= FfOrder) {
   }
   if (i != 1) {
      nextno = x * FfOrder;
   }

   // Make the seed vector
   FfSetField(basis->Field);
   FfSetNoc(basis->Noc);
   FfMulRow(vec,FF_ZERO);
   for (j = 0, x = nextno; x != 0 && j < basis->Nor; ++j, x /= FfOrder) {
      FEL co = FfFromInt(x % FfOrder);
      if (co != FF_ZERO) {
         FfAddMulRow(vec,MatGetPtr(basis,j),co);
      }
   }

   // Check for overflow
   if (x != 0) {
      return -1;
   }

   return nextno;
}
