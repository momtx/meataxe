////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Basic permutation functions
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


#define PERM_MAGIC 0x30f8326b

/// @defgroup perm Permutations
/// @details
/// In the MeatAxe, a permutation of degree n operates on {0,1,...,n-1} and is represented
/// by a Perm_t structure.
/// However, in the textual representation produced by permPrint() or by the
/// @ref prog_zpr "zpr" program, the points are numbered from 1...n.
///
/// Only permutations of equal degree can be multiplied. This can be confusing
/// because the textual representation produced by permPrint() does not include the
/// degree, and fixed points are always suppressed. For example "(3 4)(5 6 7)" could
/// be interpreted as a permutation of degree 8 or any higher degree. All these permutations
/// are in some natural way equal to each other, but they are different and incompatible
/// in the MeatAxe.
///
/// Permutations are usually created with permAlloc() or read from a file with permRead().
/// When a permutation is no longer used, the application must release the associated memory
/// by calling permFree().
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @class Perm_t
/// @details
/// Internally, a permutation is represented as an array of 32-bit integers containing the
/// images of 0,1,...,n-1. The maximum degree is 2^32-1.

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Check a permutation.
/// This function checks if the argument is a pointer to a valid permutation.
/// If the permutation is o.k., the function returns 1.
/// Otherwise, an error is signalled and, if the error handler does not
/// terminate the program, the function returns 0.
/// @return 1 if @a p is a valid permutation, 0 otherwise.

int permIsValid(const Perm_t *p)
{
   if (p == NULL) {
      return 0;
   }
   if (p->Magic != PERM_MAGIC || p->Degree < 0 || p->Data == NULL) {
      return 0;
   }
   for (int i = 0; i < p->Degree; ++i) {
      if (p->Data[i] < 0 || p->Data[i] >= p->Degree) {
         return 0;
      }
   }

   return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Checks if the given permuation is valid and aborts the program if the test fails.

void permValidate(const struct MtxSourceLocation* src, const Perm_t *p)
{
   if (p == NULL) {
      mtxAbort(src,"NULL permutation");
   }
   if (p->Magic != PERM_MAGIC || p->Degree < 0 || p->Data == NULL) {
      mtxAbort(src,"Invalid permutation (magic=%d, deg=%d)", p->Magic, p->Degree);
   }
   for (int i = 0; i < p->Degree; ++i) {
      if (p->Data[i] < 0 || p->Data[i] >= p->Degree) {
         mtxAbort(src,"Invalid value %d in permutation (deg = %d)", (int) p->Data[i], p->Degree);
      }
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Allocate a permutation
/// This function creates a permutation of the specified degree.
/// The new permutation is initialized to the identity.
/// @param deg Degree.
/// @return Pointer to the new permutation, or 0 on error.

Perm_t *permAlloc(uint32_t deg)
{
   Perm_t *p;
   int i;

   if (deg < 0) {
      mtxAbort(MTX_HERE,"deg=%d: %s",deg,MTX_ERR_BADARG);
      return NULL;
   }

   p = ALLOC(Perm_t);
   if (p == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate Perm_t structure");
      return NULL;
   }
   p->Magic = PERM_MAGIC;
   p->Degree = deg;
   p->Data = NALLOC(uint32_t,deg);
   if (p->Data == NULL) {
      sysFree(p);
      mtxAbort(MTX_HERE,"Cannot allocate permutation data");
      return NULL;
   }
   for (i = 0; i < deg; ++i) {
      p->Data[i] = i;
   }
   return p;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Free a permutation.
/// This function deletes a permutation and returns the memory to the
/// system. Do not use sysFree() on permutations because this would only
/// free the Perm_t structure but not the data buffer.
/// @param p Pointer to the permutation.
/// @return 0 on success, -1 on error.

int permFree(Perm_t *p)
{
   permValidate(MTX_HERE, p);

   sysFree(p->Data);
   memset(p,0,sizeof(Perm_t));
   sysFree(p);
   return 0;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
