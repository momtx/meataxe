////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Polymorphic objects
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

#define MAT_MAGIC 0x6233af91    /* HACK: Duplicated from  matcore.c */
#define IS_MATRIX(x) (((Matrix_t *)(x))->Magic == MAT_MAGIC)

////////////////////////////////////////////////////////////////////////////////////////////////////

static void *XRead(FILE *f)
{
   long fl;

   if (SysReadLong(f,&fl,1) != 1) {
      MTX_ERROR("Error reading file header: %S");
      return NULL;
   }
   SysFseek(f,0);
   if (fl >= 2) {
      return (void *)MatRead(f);
   } else {
      return (void *)PermRead(f);
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

   f = SysFopen(fn,FM_READ);
   if (f == NULL) {
      MTX_ERROR1("Cannot open %s: %S",fn);
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
      return MatSave((Matrix_t *)a,fn);
   }
   return PermSave((Perm_t *)a,fn);
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
      MatMul((Matrix_t *)a,(Matrix_t *)b);
   } else {
      PermMul((Perm_t *)a,(Perm_t *)b);
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
      return MatOrder((Matrix_t *)a);
   } else {
      return PermOrder((Perm_t *)a);
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
   return IS_MATRIX(a) ? (void *) MatDup(a) : (void *)PermDup(a);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Inverse of a matrix or permutation.
/// @param a The matrix or permutation.
/// @return Pointer to the inverse of a, or NULL on error.

void *XInverse(void *a)
{
   return IS_MATRIX(a) ? (void *) MatInverse(a) : (void *)PermInverse(a);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Free a matrix or permutation.
/// @param a The matrix or permutation.

void XFree(void *a)
{
   if (IS_MATRIX(a)) {
      MatFree((Matrix_t *)a);
   } else {
      PermFree((Perm_t *)a);
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
      b = (void *) MatPower((Matrix_t *) a,n);
   } else {
      b = (void *) PermPower((Perm_t *) a,n);
   }
   if (neg) {
      XFree(a);
   }
   return b;
}
