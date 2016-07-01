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

static int ErrorFlag = 0;

////////////////////////////////////////////////////////////////////////////////////////////////////

static void MyErrorHandler(const MtxErrorRecord_t *err)
{
   ErrorFlag = 1;
   err = NULL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int CheckError()
{
   int i = ErrorFlag;
   ErrorFlag = 0;
   return i;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

#define NPOLY 5

test_F FactoredPolynomialAllocation()
{
   FPoly_t *p[NPOLY];
   MtxErrorHandler_t *old_err_handler;
   int i;

   for (i = 0; i < NPOLY; ++i) {
      p[i] = FpAlloc();
   }
   for (i = 0; i < NPOLY; ++i) {
      FpIsValid(p[i]);
   }
   for (i = 0; i < NPOLY; ++i) {
      if (FpFree(p[i]) != 0) {
         Error("FpFree() failed");
      }
   }
   old_err_handler = MtxSetErrorHandler(MyErrorHandler);
   for (i = 0; i < NPOLY; ++i) {
      if (FpIsValid(p[i]) || !CheckError()) {
         Error("FpIsValid() failed");
      }
   }
   MtxSetErrorHandler(old_err_handler);
}
