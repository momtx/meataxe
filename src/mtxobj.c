////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Polymorphic objects
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


#define MAT_MAGIC 0x6233af91    /* HACK: Duplicated from  matcore.c */
#define IS_MATRIX(x) (((Matrix_t *)(x))->Magic == MAT_MAGIC)

////////////////////////////////////////////////////////////////////////////////////////////////////

static void *XRead(FILE *f)
{
   long fl;

   if (sysReadLong32(f,&fl,1) != 1) {
      mtxAbort(MTX_HERE,"Error reading file header: %S");
      return NULL;
   }
   sysFseek(f,0);
   if (fl >= 2) {
      return (void *)matRead(f);
   } else {
      return (void *)permRead(f);
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Load a matrix or permutation.
/// This function reads a matrix or permutation from a file.
/// @param fn File name.
/// @return Pointer to the matrix or permutation, NULL on error.
/// @see MatLoad PermLoad

void *XLoad(const char *fn)
{
   FILE *f;
   void *x;

   f = sysFopen(fn,"rb");
   if (f == NULL) {
      mtxAbort(MTX_HERE,"Cannot open %s: %S",fn);
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

int XSave(void *a, const char *fn)
{
   if (IS_MATRIX(a)) {
      return matSave((Matrix_t *)a,fn);
   }
   return permSave((Perm_t *)a,fn);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Multiply matrices or permutations.
/// This function multiplies a matrix or permutation by a sencond matrix
/// or permutation. Note: both objects must be of the same type.
/// @param a First matrix or permutation.
/// @param b Second matrix or permutation.

void XMul(void *a, void *b)
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

long XOrder(void *a)
{
   if (IS_MATRIX(a)) {
      return matOrder((Matrix_t *)a);
   } else {
      return permOrder((Perm_t *)a);
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Check object compatibility.
/// This function checks if two objects are compatible for |XMul()|, i.e.,
/// if they are of the same type (matrix or permutation) and have the same
/// attributes.
/// @param a First matrix or permutation.
/// @param b Second matrix or permutation.

int XIsCompatible(void *a, void *b)
{
   if (IS_MATRIX(a)) {
      Matrix_t *am = (Matrix_t *) a;
      Matrix_t *bm = (Matrix_t *) b;
      return IS_MATRIX(b) && am->Field == bm->Field && am->Nor == bm->Nor
             && am->Noc == bm->Noc;
   } else {
      Perm_t *ap = (Perm_t *) a;
      Perm_t *bp = (Perm_t *) b;
      return !IS_MATRIX(b) && ap->Degree == bp->Degree;
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Duplicate a matrix or permutation.
/// This function creates a copy of a matrix or permutation.
/// @param a The matrix or permutation.
/// @return Pointer to a copy of a, or NULL on error.
/// @see

void *XDup(void *a)
{
   return IS_MATRIX(a) ? (void *) matDup(a) : (void *)permDup(a);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Inverse of a matrix or permutation.
/// @param a The matrix or permutation.
/// @return Pointer to the inverse of a, or NULL on error.

void *XInverse(void *a)
{
   return IS_MATRIX(a) ? (void *) matInverse(a) : (void *)permInverse(a);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Free a matrix or permutation.
/// @param a The matrix or permutation.

void XFree(void *a)
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

void *XPower(void *a, int n)
{
   void *b;
   int neg = 0;

   if (n < 0) {
      a = XInverse(a);
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
      XFree(a);
   }
   return b;
}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
