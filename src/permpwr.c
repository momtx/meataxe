////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Power of a permutation
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup perm
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Power of a permutation
/// This function calculates the n-th power of a permutation.
/// It allocates a new permutation, leaving the original
/// permutation intact. The caller is responsible for deleting the
/// result when it is no longer needed.
/// @param p Pointer to the permutation.
/// @param n Exponent. Must be greather than or equal to 0.
/// @return @em n-th power of @em p or 0 on error.

Perm_t *permPower(const Perm_t *p, int n)
{

   // check arguments
   permValidate(MTX_HERE, p);
   if (n < 0) {
      mtxAbort(MTX_HERE,"Invalid exponent %d < 0",n);
      return NULL;
   }

   Perm_t* q = permAlloc(p->Degree);
   
   uint32_t* xp = p->Data;
   uint32_t* xq = q->Data;

   // calculate the n-th power
   for (uint32_t i = 0; i < p->Degree; ++i) {
      uint32_t k = i;
      for (uint32_t l = n; l > 0; --l) {
         k = xp[k];
      }
      xq[i] = k;
   }
   return q;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
