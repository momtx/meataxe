////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Test functions for matrix sets.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "check.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestMsClean1(MatrixSet_t *set)
{
   Matrix_t *zero;
   const int nor = 5;
   const int noc = 4;
   int i;

   for (i = 0; i < nor * noc; ++i) {
      int k;
      Matrix_t *m = MatAlloc(FfOrder,nor,noc);
      for (k = 0; k <= i; ++k) {
         FfInsert(MatGetPtr(m,k / noc),k % noc,FTab[k % (FfOrder - 1) + 1]);
      }
      ASSERT_EQ_INT(MsCleanAndAppend(set,m), 0);
   }

   zero = MatAlloc(FfOrder,nor,noc);
   for (i = 0; i < nor * noc; ++i) {
      Matrix_t *m = RndMat(FfOrder,nor,noc);
      MsClean(set,m);
      ASSERT_EQ_INT(MatCompare(m,zero), 0);
      MatFree(m);
   }
   MatFree(zero);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F MatrixSetClean()
{
   while (NextField() > 0) {
      MatrixSet_t *set;
      set = MsAlloc();
      TestMsClean1(set);
      MsFree(set);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F MatrixSetAllocation()
{
   while (NextField() > 0) {
      Matrix_t *m1, *m2;
      unsigned long magic;
      MatrixSet_t *set = MsAlloc();
      ASSERT(set != NULL);
      m1 = RndMat(FfOrder,10,20);
      m2 = RndMat(FfOrder,10,20);
      magic = m1->Magic;
      ASSERT(MsCleanAndAppend(set,m1) == 0);
      ASSERT(MsCleanAndAppend(set,m2) == 0);
      ASSERT(MsFree(set) == 0);
      ASSERT(m1->Magic != magic);
      ASSERT(m2->Magic != magic);
   }
}
