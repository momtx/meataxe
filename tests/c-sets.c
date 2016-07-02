////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests for integer sets
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "check.h"
#include "meataxe.h"

#include <string.h>

MTX_DEFINE_FILE_INFO

////////////////////////////////////////////////////////////////////////////////////////////////////

static int ErrorFlag = 0;

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

#define NMAT 5

test_F SetAllocaction()
{
   Set_t  *m[NMAT];
   MtxErrorHandler_t *old_err_handler;
   int i;

   for (i = 0; i < NMAT; ++i) {
      m[i] = SetAlloc();
   }

   for (i = 0; i < NMAT; ++i) {
      SetIsValid(m[i]);
      MTX_VERIFY(m[i]->Size == 0);
   }
   for (i = 0; i < NMAT; ++i) {
      if (SetFree(m[i]) != 0) {
         TST_FAIL("SetFree() failed");
      }
   }
   old_err_handler = MtxSetErrorHandler(MyErrorHandler);
   for (i = 0; i < NMAT; ++i) {
      if (SetIsValid(m[i]) || !CheckError()) {
         TST_FAIL("SetIsValid() failed");
      }
   }
   MtxSetErrorHandler(old_err_handler);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F SetBasicOperations()
{
   Set_t *s;
   long d[100];
   int i;

   memset(d,0,sizeof(d));
   for (i = 1; i <= 100; ++i) {
      int p;
      for (p = MtxRandomInt(100); d[p] != 0; p = (p + 1) % 100) {
      }
      d[p] = i;
   }

   s = SetAlloc();
   for (i = 0; i < 100; ++i) {
      int k;
      SetInsert(s,d[i]);
      if (s->Size != i + 1) {
         TST_FAIL("Bad size");
      }
      for (k = 0; k <= i; ++k) {
         if (!SetContains(s,d[k])) {
            TST_FAIL("Element not inserted");
         }
      }
      for (k = i + 1; k < 100; ++k) {
         if (SetContains(s,d[k])) {
            TST_FAIL("Unexpected Element");
         }
      }
   }
   SetFree(s);
}
