////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check functions for factored polynomials.
//
// (C) Copyright 1998-2016 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "check.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

#define NPOLY 5

test_F FactoredPolynomialAllocation()
{
   FPoly_t *p[NPOLY];
   int i;

   for (i = 0; i < NPOLY; ++i) {
      p[i] = FpAlloc();
   }
   for (i = 0; i < NPOLY; ++i) {
      FpIsValid(p[i]);
   }
   for (i = 0; i < NPOLY; ++i) {
      ASSERT(FpFree(p[i]) == 0);
   }
   TstStartErrorChecking();
   for (i = 0; i < NPOLY; ++i) {
      ASSERT(!FpIsValid(p[i]) && TstHasError());
   }
}
