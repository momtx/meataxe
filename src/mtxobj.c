////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Polymorphic objects
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>
#include <errno.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


#define IS_MATRIX(x) (((Matrix_t *)x)->typeId == MTX_TYPE_MATRIX)

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Load a matrix or permutation.
/// This function reads a matrix or permutation from a file.
/// @param fn File name.
/// @return Pointer to the matrix or permutation, NULL on error.
/// @see MatLoad PermLoad

void* objLoad(const char* fn)
{
   MtxFile_t* f = mfOpen(fn, "rb");
   const uint32_t objectType = mfReadHeader(f);
   void* x = NULL;
   switch (objectType) {
      case MTX_TYPE_MATRIX:
         x = matReadData(f);
         break;
      case MTX_TYPE_PERMUTATION:
         x = permReadData(f);
         break;
      default:
         mtxAbort(MTX_HERE, "%s: unsuported object type 0x%lx",
            f->name, (unsigned long)objectType);
   }

   mfClose(f);
   return x;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Save a matrix or permutation.
/// This function writes a matrix or permutation to a file. If a file with
/// the given name exists, it is destroyed.
/// @param a Matrix or permutation.
/// @param fn File name.
/// @return $0$ on success, $-1$ on error.

int objSave(void *a, const char *fn)
{
   if (IS_MATRIX(a)) {
      matSave((Matrix_t *)a,fn);
   } else {
      permSave((Perm_t *)a,fn);
   }
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Multiply matrices or permutations.
/// This function multiplies a matrix or permutation by a sencond matrix
/// or permutation. Note: both objects must be of the same type.
/// @param a First matrix or permutation.
/// @param b Second matrix or permutation.

void objMul(void *a, void *b)
{
   if (IS_MATRIX(a)) {
      matMul((Matrix_t *)a,(Matrix_t *)b);
   } else {
      permMul((Perm_t *)a,(Perm_t *)b);
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns the order of a matrix or permutation.

long objOrder(void* a)
{
   if (IS_MATRIX(a)) {
      return matOrder((Matrix_t*)a);
   }
   else {
      return permOrder((Perm_t*)a);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Check object compatibility.
///
/// This function checks if two objects can be multiplied by @ref objMul.

int objCanMultiply(void* a, void* b)
{
   if (IS_MATRIX(a)) {
      Matrix_t* am = (Matrix_t*) a;
      Matrix_t* bm = (Matrix_t*) b;
      return IS_MATRIX(b) && am->field == bm->field && am->nor == bm->nor
             && am->noc == bm->noc;
   }
   else {
      Perm_t* ap = (Perm_t*) a;
      Perm_t* bp = (Perm_t*) b;
      return !IS_MATRIX(b) && ap->degree == bp->degree;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Create an independent copy of a matrix or permutation.

void* objDup(void* a)
{
   return IS_MATRIX(a) ? (void*)matDup(a) : (void*)permDup(a);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns the inverse of a matrix or permutation.

void *objInverse(void *a)
{
   return IS_MATRIX(a) ? (void *) matInverse(a) : (void *)permInverse(a);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Destroys a matrix or permutation.

void objFree(void *a)
{
   if (IS_MATRIX(a)) {
      matFree((Matrix_t *)a);
   } else {
      permFree((Perm_t *)a);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Calculates the n-th power of a matrix or permutation. @a n may be negative if the object is
/// invertible.

void *objPower(void *a, int n)
{
   void *b;
   int neg = 0;

   if (n < 0) {
      a = objInverse(a);
      if (a == NULL) {
         return NULL;
      }
      n = -n;
      neg = 1;
   }

   if (IS_MATRIX(a)) {
      b = (void *) matPower((Matrix_t *) a,n);
   } else {
      b = (void *) permPower((Perm_t *) a,n);
   }
   if (neg) {
      objFree(a);
   }
   return b;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
