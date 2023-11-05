////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Map a vector under a permutation.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

/// @addtogroup ff
/// @{
 
/// Multiplies a vector from the right by a permutation.
///
/// This function multiplies the vector @a row from the right with the permutation @a perm
/// and stores the result in @a result. More explicitly: if <tt>perm[i] = k</tt>, then the
/// i-th mark of the vector is stored in the k-th position of the result. 
///
/// @param result Result vector (@a noc columns).
/// @param row A row vector with @a noc columns.
/// @param perm Pointer to a table of @a noc numbers defining a permutation of {0,..., noc-1}.
/// @param noc Number of columns in @a row and @result
///
/// Note: @a result and @a row must not overlap.

void ffPermRow(PTR result, PTR row, const uint32_t* perm, int noc)
{
   register FEL f;
   register int i;

   // Verify that row and result do not overlap
   MTX_ASSERT(row != result);
   MTX_ASSERT_DEBUG(ffGetPtr(row, 1, noc) <= result || row >= ffGetPtr(result,1, noc));

   const uint32_t *p = perm;
   for (i = 0; i < noc; ++i)
   {
      MTX_ASSERT_DEBUG(*p >= 0 && *p < noc);
      f = ffExtract(row,i);
      ffInsert(result,*p++,f);
   }
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
