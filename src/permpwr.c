////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Power of a permutation
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

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

Perm_t *PermPower(const Perm_t *p, int n)
{
   Perm_t *q;
   long *xp;
   long *xq;
   long i, k, l;

   // Check arguments
   if (!PermIsValid(p)) {
      return NULL;
   }
   if (n < 0) {
      MTX_ERROR1("Invalid exponent %d < 0",n);
      return NULL;
   }

   // Allocate a new permutation for the result
   q = PermAlloc(p->Degree);
   if (q == NULL) {
      return NULL;
   }
   xp = p->Data;
   xq = q->Data;

   // Calculate the n-th power
   for (i = 0; i < p->Degree; ++i) {
      for (k = i, l = n; l > 0; --l) {
         k = xp[k];
      }
      xq[i] = k;
   }
   return q;
}


/// @}