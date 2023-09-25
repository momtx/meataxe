////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Polymorphic objects
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>
#include <errno.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


#define IS_MATRIX(x) (((Matrix_t *)x)->typeId == MTX_TYPE_MATRIX)
#define IS_PERMUTATION(x) (((Perm_t *)x)->typeId == MTX_TYPE_PERMUTATION)

////////////////////////////////////////////////////////////////////////////////////////////////////

static void *XRead(FILE *f)
{
   uint32_t header[3];
   sysRead32(f,header,3);
   if (header[0] < MTX_TYPE_BEGIN)
      return matReadData(f, header);
   if (header[0] == MTX_TYPE_PERMUTATION)
      return permReadData(f, header);
   mtxAbort(MTX_HERE, "Unsupported object type 0x%lx", (unsigned long) header[0]);
   return NULL;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Load a matrix or permutation.
/// This function reads a matrix or permutation from a file.
/// @param fn File name.
/// @return Pointer to the matrix or permutation, NULL on error.
/// @see MatLoad PermLoad

void *objLoad(const char *fn)
{
   FILE *f;
   void *x;

   f = sysFopen(fn,"rb");
   if (f == NULL) {
      mtxAbort(MTX_HERE,"Cannot open %s: %s",fn, strerror(errno));
      return NULL;
   }
   x = XRead(f);
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
/// Order of a matrix or permutation.
/// This function calculates the order of a matrix or permutation.
/// @param a The matrix or permutation.
/// @see MatOrder PermOrder

long objOrder(void *a)
{
   if (IS_MATRIX(a)) {
      return matOrder((Matrix_t *)a);
   } else {
      return permOrder((Perm_t *)a);
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Check object compatibility.
/// This function checks if two objects are compatible for |objMul()|, i.e.,
/// if they are of the same type (matrix or permutation) and have the same
/// attributes.
/// @param a First matrix or permutation.
/// @param b Second matrix or permutation.

int objCanMultiply(void *a, void *b)
{
   if (IS_MATRIX(a)) {
      Matrix_t *am = (Matrix_t *) a;
      Matrix_t *bm = (Matrix_t *) b;
      return IS_MATRIX(b) && am->field == bm->field && am->nor == bm->nor
             && am->noc == bm->noc;
   } else {
      Perm_t *ap = (Perm_t *) a;
      Perm_t *bp = (Perm_t *) b;
      return !IS_MATRIX(b) && ap->degree == bp->degree;
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Duplicate a matrix or permutation.
/// This function creates a copy of a matrix or permutation.
/// @param a The matrix or permutation.
/// @return Pointer to a copy of a, or NULL on error.
/// @see

void *objDup(void *a)
{
   return IS_MATRIX(a) ? (void *) matDup(a) : (void *)permDup(a);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Inverse of a matrix or permutation.
/// @param a The matrix or permutation.
/// @return Pointer to the inverse of a, or NULL on error.

void *objInverse(void *a)
{
   return IS_MATRIX(a) ? (void *) matInverse(a) : (void *)permInverse(a);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Free a matrix or permutation.
/// @param a The matrix or permutation.

void objFree(void *a)
{
   if (IS_MATRIX(a)) {
      matFree((Matrix_t *)a);
   } else {
      permFree((Perm_t *)a);
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Power of a matrix or permutation.
/// This function calculates the n-th power of a matrix or permutation.
/// $n$ may be negative.
/// @param a The matrix or permutation.
/// @param n The power.
/// @return n-th power of a, or NULL on error.

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
